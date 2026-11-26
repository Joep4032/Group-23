function results = faRunLoadStudy(Cyl, o, mainTabs)
%FARUNLOADSTUDY Fuel-mass adjustment at prescribed gross IMEP.
%
% Inputs:
%   Cyl      Existing engine geometry.
%   o        Options returned by faOptions().
%   mainTabs Existing tab-group handle, or [] for calculations only.
%
% Requires:
%   faSimulate.m
%   faFuels.m
%   faMixture.m
%   faWiebe.m
%   Existing CylinderVolume and NASA functions.
%
% All new plots appear inside mainTabs.
% No additional figure window is created.

if nargin < 3
    mainTabs = [];
end

%% Locate project folders

projectRoot = fileparts(fileparts(mfilename('fullpath')));

originalPath = path;
pathCleanup = onCleanup(@() path(originalPath)); %#ok<NASGU>

addpath(fullfile(projectRoot, 'Functions'));
addpath(fullfile(projectRoot, 'Nasa'));

%% Load-study settings

targetLoads = [1.5 3.0 4.5];       % Gross IMEP [bar]
loadTolerance = 5e-4;             % Acceptable error [bar]

assert(o.RPM == 1500, ...
    'The load study requires 1500 RPM.');

assert(all(o.SOI_bTDC_deg >= 14 & o.SOI_bTDC_deg <= 21), ...
    'Injection timing must remain between 14 and 21 degrees BTDC.');

assert(all(ismember(diff(o.SOI_bTDC_deg), [2 3])), ...
    'SOI values must increase in steps of 2 or 3 degrees.');

if ~isempty(mainTabs)
    assert(isgraphics(mainTabs, 'uitabgroup'), ...
        'mainTabs must be an existing tab-group handle.');
end

%% Load NASA species

global Runiv

if isempty(Runiv)
    Runiv = 8.314;
end

database = load( ...
    fullfile(projectRoot, 'Nasa', 'NasaThermalDatabase.mat'), ...
    'Sp');

speciesNames = {'O2', 'N2', 'CO2', 'H2O'};
indices = zeros(1, 4);

for k = 1:4
    index = find(strcmp({database.Sp.Name}, speciesNames{k}), 1);

    assert(~isempty(index), ...
        'NASA species not found: %s', speciesNames{k});

    indices(k) = index;
end

species = database.Sp(indices);

%% Reuse existing fuels and engine geometry

fuels = faFuels();

nFuel = numel(fuels);
nLoad = numel(targetLoads);
nTiming = numel(o.SOI_bTDC_deg);

displacement = pi*Cyl.Bore^2*Cyl.Stroke/4;
cyclesPerSecond = o.RPM/120;

% BDC volume, accounting for the geometry function's TDC offset.
V_BDC = CylinderVolume(Cyl.TDCangle-180, Cyl);

% Initial oxygen availability limits the admissible fuel mass.
nAir = o.P_BDC_Pa*V_BDC/(Runiv*o.T_BDC_K);
nO2available = nAir*o.AirO2MoleFraction;

% Atomic masses consistent with the supplied species database.
M = [species.Mass];

MO = M(1)/2;
MH = (M(4)-MO)/2;
MC = M(3)-2*MO;

%% Storage

cycles = cell(nFuel, nLoad, nTiming);
rows = cell(nFuel*nLoad*nTiming, 1);

rowIndex = 0;

fprintf('\nStarting load-controlled theoretical study.\n');
fprintf('Targets: 1.5, 3.0 and 4.5 bar GROSS IMEP.\n');
fprintf('Fuel mass is adjusted at every injection timing.\n');

%% Fuel, load and timing sweep

for iFuel = 1:nFuel

    fuel = fuels(iFuel);

    carbon = fuel.CHO(1);
    hydrogen = fuel.CHO(2);
    oxygen = fuel.CHO(3);

    fuelMolarMass = carbon*MC + hydrogen*MH + oxygen*MO;
    oxygenDemand = carbon + hydrogen/4 - oxygen/2;

    assert(oxygenDemand > 0 && fuelMolarMass > 0, ...
        'Invalid elemental composition for %s.', fuel.Name);

    % Maximum fuel mass with a small excess-oxygen margin.
    maximumMass = ...
        0.98*nO2available/oxygenDemand*fuelMolarMass;

    minimumMass = 1e-9;

    assert(maximumMass > minimumMass, ...
        'No admissible fuel-mass interval for %s.', fuel.Name);

    for iLoad = 1:nLoad

        targetIMEP = targetLoads(iLoad);

        for iTiming = 1:nTiming

            soi = o.SOI_bTDC_deg(iTiming);

            %% Define the load residual

            objective = @(mass) loadResidual( ...
                mass, fuel, soi, targetIMEP);

            lowerMass = minimumMass;

            upperMass = min( ...
                maximumMass, ...
                max(o.ReferenceFuelMass_kg, 1e-6));

            lowerResidual = objective(lowerMass);
            upperResidual = objective(upperMass);

            % Expand the bracket without exceeding oxygen availability.
            while upperResidual < 0 && upperMass < maximumMass

                upperMass = min(1.5*upperMass, maximumMass);
                upperResidual = objective(upperMass);

            end

            assert(lowerResidual < 0 && upperResidual >= 0, ...
                ['Cannot reach %.1f bar gross IMEP for %s at ', ...
                 '%g degrees BTDC within the lean-combustion range.'], ...
                targetIMEP, fuel.Name, soi);

            %% Solve for the required fuel mass

            rootOptions = optimset( ...
                'Display', 'off', ...
                'TolX', 1e-12);

            [fuelMass, ~, exitFlag] = fzero( ...
                objective, ...
                [lowerMass upperMass], ...
                rootOptions);

            assert(exitFlag > 0, ...
                'Fuel-mass adjustment did not converge.');

            %% Run the final cycle at the solved mass

            caseOptions = o;

            % This makes faSimulate use the mass supplied below directly.
            caseOptions.ComparisonMode = 'equalMass';
            caseOptions.ReferenceFuelMass_kg = fuelMass;

            c = faSimulate( ...
                Cyl, fuel, soi, caseOptions, species);

            achievedIMEP = c.Work_J/displacement/1e5;
            loadError = achievedIMEP-targetIMEP;

            assert(abs(loadError) <= loadTolerance, ...
                ['Load matching failed for %s at %g BTDC: ', ...
                 'target %.4f bar, achieved %.4f bar.'], ...
                fuel.Name, soi, targetIMEP, achievedIMEP);

            cycles{iFuel, iLoad, iTiming} = c;

            %% KPI calculations

            fuelFlow = c.FuelMass_kg*cyclesPerSecond;
            co2Flow = c.CO2Mass_kg*cyclesPerSecond;
            indicatedPower = c.Work_J*cyclesPerSecond;

            row = struct();

            row.Fuel = string(fuel.Name);
            row.RPM = o.RPM;

            row.TargetGrossIMEP_bar = targetIMEP;
            row.GrossIMEP_bar = achievedIMEP;
            row.IMEPerror_bar = loadError;

            row.SOI_bTDC_deg = soi;

            row.GrossIndicatedEfficiency_pct = ...
                100*c.Work_J/c.FuelEnergy_J;

            row.BrakeEfficiency_pct = NaN;

            row.GrossISFC_gkWh = ...
                c.FuelMass_kg/c.Work_J*3.6e9;

            row.BSFC_gkWh = NaN;

            row.GrossISCO2_gkWh = ...
                c.CO2Mass_kg/c.Work_J*3.6e9;

            row.BSCO2_gkWh = NaN;

            row.Fuel_mgCycle = c.FuelMass_kg*1e6;
            row.FuelEnergy_JCycle = c.FuelEnergy_J;
            row.FuelFlow_gs = fuelFlow*1000;

            row.GrossIndicatedPower_kW = indicatedPower/1000;
            row.BrakePower_kW = NaN;

            row.GrossWork_JCycle = c.Work_J;
            row.PeakPressure_bar = max(c.Pressure_Pa)/1e5;

            row.CA50_degATDC = c.CA50_degATDC;
            row.SOC_degATDC = c.SOC_degATDC;
            row.IgnitionDelay_ms = c.IgnitionDelay_ms;

            row.Lambda = c.Lambda;
            row.ExhaustCO2Flow_gs = co2Flow*1000;

            row.EnergyBalanceError_pct = ...
                c.EnergyBalanceError_pct;

            rowIndex = rowIndex+1;
            rows{rowIndex} = row;

            fprintf( ...
                '%2d/%2d | %-9s | %.1f bar | %2g BTDC | %.3f mg/cycle\n', ...
                rowIndex, numel(rows), fuel.Name, ...
                targetIMEP, soi, c.FuelMass_kg*1e6);

        end
    end
end

%% Assemble outputs

results = struct();

results.KPI = struct2table(vertcat(rows{:}));
results.Cycles = cycles;
results.TargetGrossIMEP_bar = targetLoads;
results.Fuels = fuels;
results.Options = o;
results.Geometry = Cyl;
results.LoadTolerance_bar = loadTolerance;
results.OutputDirectory = '';

%% Export results

if o.ExportResults

    outputRoot = fullfile(projectRoot, 'LoadStudyResults');

    if ~exist(outputRoot, 'dir')
        mkdir(outputRoot);
    end

    % Preserve earlier studies by creating a unique output directory.
    outputDirectory = tempname(outputRoot);
    mkdir(outputDirectory);

    results.OutputDirectory = outputDirectory;

    writetable(results.KPI, ...
        fullfile(outputDirectory, 'All_Load_KPIs.csv'));

    for iLoad = 1:nLoad

        selected = abs( ...
            results.KPI.TargetGrossIMEP_bar-targetLoads(iLoad)) < 1e-9;

        loadTable = results.KPI(selected,:);

        filename = sprintf( ...
            'Load_%.1fbar_KPIs.csv', targetLoads(iLoad));

        writetable(loadTable, ...
            fullfile(outputDirectory, filename));

    end

    save(fullfile(outputDirectory, 'LoadStudyResults.mat'), ...
        'results');

end

%% Plot in the existing window

if o.MakePlots && ~isempty(mainTabs)

    % Remove only pages created by a previous run of this function.
    delete(findall(mainTabs, ...
        'Type', 'uitab', ...
        'Tag', 'FuelLoadStudy'));

    fuelColors = lines(nFuel);
    timingColors = lines(nTiming);

    fuelLabels = {fuels.Name};

    timingLabels = arrayfun( ...
        @(x) sprintf('%g deg BTDC', x), ...
        o.SOI_bTDC_deg, ...
        'UniformOutput', false);

    % Use a common pressure scale across all loads and fuels.
    pressureMaximum = max(results.KPI.PeakPressure_bar);
    pressureLimit = 10*ceil(pressureMaximum/10);

    volumeLimit = ceil(V_BDC*1000*10)/10;

    for iLoad = 1:nLoad

        targetIMEP = targetLoads(iLoad);

        selected = abs( ...
            results.KPI.TargetGrossIMEP_bar-targetIMEP) < 1e-9;

        loadTable = results.KPI(selected,:);

        %% One parent page per load, containing several subpages

        loadPage = uitab(mainTabs, ...
            'Title', sprintf('Load %.1f bar', targetIMEP), ...
            'Tag', 'FuelLoadStudy');

        loadTabs = uitabgroup('Parent', loadPage);

        %% pV timing sweep

        pvPage = uitab(loadTabs, 'Title', 'SOI pV diagrams');

        for iFuel = 1:nFuel

            % CHANGED 25/11: grid grows with the number of fuels (was 2x2,
            % which fails with the two added B7-OME3 blends: 6 fuels).
            ax = subplot(2, ceil(nFuel/2), iFuel, 'Parent', pvPage);
            hold(ax, 'on');

            for iTiming = 1:nTiming

                c = cycles{iFuel, iLoad, iTiming};

                % Idealized constant-volume closure; no pumping loop.
                volume = [c.Volume_m3; c.Volume_m3(1)]*1000;
                pressure = [c.Pressure_Pa; c.Pressure_Pa(1)]/1e5;

                plot(ax, volume, pressure, ...
                    'Color', timingColors(iTiming,:), ...
                    'LineWidth', 1.4);

            end

            xlabel(ax, 'Volume [liter]');
            ylabel(ax, 'Pressure [bar absolute]');

            title(ax, sprintf('%s: %.1f bar gross IMEP', ...
                fuels(iFuel).Name, targetIMEP));

            xlim(ax, [0 volumeLimit]);
            ylim(ax, [0 pressureLimit]);

            grid(ax, 'on');
            legend(ax, timingLabels, 'Location', 'best');

        end

        %% Pressure and heat release at the reference timing

        referenceIndex = find( ...
            abs(o.SOI_bTDC_deg-o.ReferenceSOI_bTDC_deg) < 1e-9, 1);

        assert(~isempty(referenceIndex), ...
            'Reference timing must be included in the SOI sweep.');

        heatPage = uitab(loadTabs, ...
            'Title', 'Pressure and heat release');

        axPressure = subplot(2,1,1, 'Parent', heatPage);
        axHeat = subplot(2,1,2, 'Parent', heatPage);

        hold(axPressure, 'on');
        hold(axHeat, 'on');

        for iFuel = 1:nFuel

            c = cycles{iFuel, iLoad, referenceIndex};

            plot(axPressure, ...
                c.Theta_degATDC, c.Pressure_Pa/1e5, ...
                'Color', fuelColors(iFuel,:), ...
                'LineWidth', 1.5);

            plot(axHeat, ...
                c.Theta_degATDC, c.HeatRelease_Jdeg, ...
                'Color', fuelColors(iFuel,:), ...
                'LineWidth', 1.5);

        end

        xlabel(axPressure, 'Crank angle relative to TDC [deg]');
        ylabel(axPressure, 'Pressure [bar absolute]');

        title(axPressure, sprintf( ...
            '%.1f bar gross IMEP; SOI %g deg BTDC', ...
            targetIMEP, o.ReferenceSOI_bTDC_deg));

        xlabel(axHeat, 'Crank angle relative to TDC [deg]');
        ylabel(axHeat, 'Prescribed heat release [J/deg]');

        xlim(axPressure, [-40 100]);
        xlim(axHeat, [-40 100]);

        grid(axPressure, 'on');
        grid(axHeat, 'on');

        legend(axPressure, fuelLabels, 'Location', 'best');
        legend(axHeat, fuelLabels, 'Location', 'best');

        %% KPI trends

        trendPage = uitab(loadTabs, 'Title', 'Timing KPIs');

        fields = {
            'GrossIndicatedEfficiency_pct'
            'Fuel_mgCycle'
            'GrossISFC_gkWh'
            'PeakPressure_bar'
            'GrossISCO2_gkWh'
            'GrossIndicatedPower_kW'
        };

        labels = {
            'Gross indicated efficiency [%]'
            'Fuel mass [mg/cycle]'
            'Gross ISFC [g/kWh]'
            'Peak pressure [bar absolute]'
            'Gross ISCO2 [g/kWh]'
            'Gross indicated power [kW]'
        };

        for iPlot = 1:numel(fields)

            ax = subplot(3,2,iPlot, 'Parent', trendPage);
            hold(ax, 'on');

            for iFuel = 1:nFuel

                fuelRows = ...
                    loadTable.Fuel == string(fuels(iFuel).Name);

                fuelTable = sortrows( ...
                    loadTable(fuelRows,:), 'SOI_bTDC_deg');

                plot(ax, ...
                    fuelTable.SOI_bTDC_deg, ...
                    fuelTable.(fields{iPlot}), ...
                    '-o', ...
                    'Color', fuelColors(iFuel,:), ...
                    'LineWidth', 1.4);

            end

            xlabel(ax, 'SOI [degrees BTDC]');
            ylabel(ax, labels{iPlot});

            grid(ax, 'on');

            if iPlot == 1
                legend(ax, fuelLabels, 'Location', 'best');
            end

        end

        %% KPI table

        tablePage = uitab(loadTabs, 'Title', 'KPI table');

        tableData = table2cell(loadTable);
        tableData(:,1) = cellstr(loadTable.Fuel);

        uitable( ...
            'Parent', tablePage, ...
            'Units', 'normalized', ...
            'Position', [0.01 0.06 0.98 0.93], ...
            'Data', tableData, ...
            'ColumnName', loadTable.Properties.VariableNames, ...
            'RowName', []);

        uicontrol( ...
            'Parent', tablePage, ...
            'Style', 'text', ...
            'Units', 'normalized', ...
            'Position', [0.01 0.005 0.98 0.045], ...
            'HorizontalAlignment', 'left', ...
            'String', ...
            ['Gross indicated load; gas exchange omitted. ', ...
             'Brake results unavailable. CO2 is not lifecycle CO2.']);

    end
end

fprintf('\nLoad study complete: %d operating points.\n', rowIndex);

fprintf('Maximum load error: %.6f bar.\n', ...
    max(abs(results.KPI.IMEPerror_bar)));

%% Nested objective used by fzero

    function residual = loadResidual(mass, fuel, soi, targetIMEP)

        trialOptions = o;

        trialOptions.ComparisonMode = 'equalMass';
        trialOptions.ReferenceFuelMass_kg = mass;

        trialCycle = faSimulate( ...
            Cyl, fuel, soi, trialOptions, species);

        predictedIMEP = trialCycle.Work_J/displacement/1e5;

        residual = predictedIMEP-targetIMEP;

    end

end