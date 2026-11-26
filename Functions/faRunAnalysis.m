function results = faRunAnalysis(Cyl, o)
%FARUNANALYSIS Run theoretical fuel and injection-timing comparisons.
% Reuses the existing engine geometry and NASA functions.
% No target IMEP is imposed. Brake-specific quantities remain unavailable.

root = fileparts(fileparts(mfilename('fullpath')));

originalPath = path;
pathCleanup = onCleanup(@() path(originalPath)); %#ok<NASGU>

addpath(fullfile(root, 'Functions'));
addpath(fullfile(root, 'Nasa'));

%% Validate engine geometry

required = {'Bore', 'Stroke', 'CompressionRatio', ...
            'ConRod', 'TDCangle'};

for k = 1:numel(required)
    assert(isfield(Cyl, required{k}), ...
        'Missing engine parameter: Cyl.%s', required{k});

    validateattributes(Cyl.(required{k}), {'numeric'}, ...
        {'scalar', 'real', 'finite'});
end

assert(Cyl.Bore > 0 && Cyl.Stroke > 0, ...
    'Bore and stroke must be positive.');

assert(Cyl.CompressionRatio > 1 && ...
       Cyl.ConRod > Cyl.Stroke/2, ...
    'Invalid compression ratio or connecting-rod length.');

%% Validate operating conditions

assert(isscalar(o.RPM) && o.RPM == 1500, ...
    'This study requires 1500 RPM.');

validateattributes(o.SOI_bTDC_deg, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonempty'});

assert(all(o.SOI_bTDC_deg >= 14 & o.SOI_bTDC_deg <= 21), ...
    'SOI must lie between 14 and 21 degrees BTDC.');

increments = diff(o.SOI_bTDC_deg);

assert(all(increments >= 2 & increments <= 3), ...
    'Use increasing SOI values with steps of 2 to 3 CAD.');

assert(any(abs(o.SOI_bTDC_deg - ...
    o.ReferenceSOI_bTDC_deg) < 1e-9), ...
    'The reference timing must be included in the SOI sweep.');

assert(any(strcmp(o.ComparisonMode, ...
    {'equalEnergy', 'equalMass'})), ...
    'ComparisonMode must be equalEnergy or equalMass.');

%% Validate model parameters

positiveFields = {
    'ReferenceFuelMass_kg', 'ReferenceLHV_Jkg', ...
    'P_BDC_Pa', 'T_BDC_K', 'ReferenceIgnitionDelay_s', ...
    'ReferenceCetaneNumber', 'BurnDuration_deg', 'WiebeA', ...
    'ReferenceTemperature_K', 'Step_deg', 'EnergyTolerance_J'
};

for k = 1:numel(positiveFields)
    validateattributes(o.(positiveFields{k}), {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
end

validateattributes(o.WiebeM, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'});

validateattributes(o.CetaneDelayExponent, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'});

validateattributes(o.UnburnedFuelHeatCapacity_JkgK, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'});

validateattributes(o.MaximumIterations, {'numeric'}, ...
    {'scalar', 'integer', 'positive'});

validateattributes(o.TemperatureBounds_K, {'numeric'}, ...
    {'vector', 'numel', 2, 'real', 'finite', 'positive'});

assert(diff(o.TemperatureBounds_K) > 0, ...
    'Temperature bounds must increase.');

assert(o.T_BDC_K >= o.TemperatureBounds_K(1) && ...
       o.T_BDC_K <= o.TemperatureBounds_K(2), ...
    'Initial temperature is outside the configured interval.');

assert(o.AirO2MoleFraction > 0 && o.AirO2MoleFraction < 1, ...
    'Invalid oxygen mole fraction.');

assert(o.Step_deg <= 1 && ...
    abs(360/o.Step_deg - round(360/o.Step_deg)) < 1e-8, ...
    'Step_deg must divide 360 exactly and be no larger than 1 CAD.');

%% Load the existing NASA database

global Runiv

if isempty(Runiv)
    Runiv = 8.314;
end

validateattributes(Runiv, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});

databasePath = fullfile(root, 'Nasa', 'NasaThermalDatabase.mat');

assert(isfile(databasePath), ...
    'NASA database not found: %s', databasePath);

database = load(databasePath, 'Sp');

speciesNames = {'O2', 'N2', 'CO2', 'H2O'};
indices = zeros(1, 4);

for k = 1:4
    index = find(strcmp({database.Sp.Name}, speciesNames{k}), 1);

    assert(~isempty(index), ...
        'Missing NASA species: %s', speciesNames{k});

    indices(k) = index;
end

species = database.Sp(indices);
fuels = faFuels();

%% Prepare the fuel and timing sweep

nFuel = numel(fuels);
nTiming = numel(o.SOI_bTDC_deg);

cycles = cell(nFuel, nTiming);

% Store individual result structures in cells.
% This avoids assigning into a structure array with incompatible fields.
rows = cell(nFuel*nTiming, 1);

displacement = pi*Cyl.Bore^2*Cyl.Stroke/4;

% Single-cylinder, four-stroke engine:
cyclesPerSecond = o.RPM/120;

rowIndex = 0;

%% Run each fuel at each injection timing

for i = 1:nFuel
    for j = 1:nTiming

        c = faSimulate( ...
            Cyl, fuels(i), o.SOI_bTDC_deg(j), o, species);

        cycles{i, j} = c;

        rowIndex = rowIndex + 1;

        % Start a fresh result structure for every simulation.
        row = struct();

        fuelFlow = c.FuelMass_kg*cyclesPerSecond;
        co2Flow = c.CO2Mass_kg*cyclesPerSecond;
        indicatedPower = c.Work_J*cyclesPerSecond;

        %% Identification

        row.Fuel = string(fuels(i).Name);
        row.RPM = o.RPM;
        row.SOI_bTDC_deg = o.SOI_bTDC_deg(j);

        %% Efficiency

        row.GrossIndicatedEfficiency_pct = ...
            100*c.Work_J/c.FuelEnergy_J;

        row.BrakeEfficiency_pct = NaN;

        %% Specific fuel consumption and CO2

        if c.Work_J > 0
            row.GrossISFC_gkWh = ...
                c.FuelMass_kg/c.Work_J*3.6e9;

            row.GrossISCO2_gkWh = ...
                c.CO2Mass_kg/c.Work_J*3.6e9;
        else
            row.GrossISFC_gkWh = NaN;
            row.GrossISCO2_gkWh = NaN;
        end

        row.BSFC_gkWh = NaN;
        row.BSCO2_gkWh = NaN;

        %% Fuel input

        row.Fuel_mgCycle = c.FuelMass_kg*1e6;
        row.FuelEnergy_JCycle = c.FuelEnergy_J;
        row.FuelFlow_gs = fuelFlow*1000;

        %% Power, work and resulting load

        row.GrossIndicatedPower_kW = indicatedPower/1000;
        row.BrakePower_kW = NaN;

        row.GrossWork_JCycle = c.Work_J;
        row.GrossIMEP_bar = c.Work_J/displacement/1e5;

        %% Pressure and combustion phasing

        row.PeakPressure_bar = max(c.Pressure_Pa)/1e5;
        row.CA50_degATDC = c.CA50_degATDC;
        row.Lambda = c.Lambda;

        %% CO2 and numerical diagnostics

        row.ExhaustCO2Flow_gs = co2Flow*1000;
        row.EnergyBalanceError_pct = c.EnergyBalanceError_pct;
        row.SOC_degATDC = c.SOC_degATDC;
        row.IgnitionDelay_ms = c.IgnitionDelay_ms;

        % Cell assignment instead of structure-array assignment.
        rows{rowIndex} = row;

    end
end

%% Assemble the results

results = struct();

% Combine the consistently constructed rows and convert to a table.
results.KPI = struct2table(vertcat(rows{:}));

results.Cycles = cycles;
results.Fuels = fuels;
results.Options = o;
results.Geometry = Cyl;
results.OutputDirectory = '';

%% Export results

if o.ExportResults

    outputRoot = fullfile(root, 'FuelAnalysisResults');

    if ~exist(outputRoot, 'dir')
        mkdir(outputRoot);
    end

    % Unique directory prevents overwriting previous results.
    outputDirectory = tempname(outputRoot);
    mkdir(outputDirectory);

    results.OutputDirectory = outputDirectory;

    writetable(results.KPI, ...
        fullfile(outputDirectory, 'FuelAnalysis_KPI.csv'));

    save(fullfile(outputDirectory, ...
        'FuelAnalysis_results.mat'), 'results');

end

%% Display the tabbed figures

if o.MakePlots
    faPlotAnalysis(results);
end

end