%% Experimental results analysis (ADDED 25/11)
% Analyses every measured run listed in Data/ExperimentRuns.csv:
% different fuels (incl. the B7-OME3 blends), loads, injection timings,
% engines and repeats, each with its own SOI, intake conditions and flows.
%
% How to use:
%   1. Put the lab pressure files (text: crank angle, pressure) in Data/.
%   2. Add one row per run to Data/ExperimentRuns.csv (repeats = own row).
%      Fuel names: Diesel (or B7), HVO, OME3, FT diesel,
%                  B7-OME3 12.5%, B7-OME3 25%   (see faFuels.m)
%      Empty cells are allowed; the KPIs that need them become NaN.
%   3. Run this script. It works on its own; Simple.m is not needed.
%
% Output: KPI table in the workspace (ExperimentKPI), a tabbed figure,
% and ExperimentResults/<date_time>/ with the CSV and .mat file.

%% Paths and NASA data
experimentRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(experimentRoot, 'Functions'), fullfile(experimentRoot, 'Nasa'));

global Runiv
Runiv = 8.314;

% Same species and order as Simple.m: SpS(2:5) = O2, N2, CO2, H2O
[SpS, ~] = myload(fullfile(experimentRoot, 'Nasa', 'NasaThermalDatabase.mat'), ...
    {'Diesel', 'O2', 'N2', 'CO2', 'H2O'});

%% Engine geometry (same values as Simple.m)
ExpCyl.Bore             = 104e-3;
ExpCyl.Stroke           = 85e-3;
ExpCyl.CompressionRatio = 21.5;
ExpCyl.ConRod           = 136.5e-3;
ExpCyl.TDCangle         = 0;

%% Options (change here if the lab files differ)
eo = struct();

% Data file layout
eo.CrankAngleColumn = 1;          % column with crank angle [deg]
eo.PressureColumn = 2;            % column with cylinder pressure
eo.PressureUnit_Pa = 1e5;         % 1e5 if the file is in bar
eo.CycleStart_deg = -360;         % first crank angle of a cycle

% Engine events and defaults for missing inputs
eo.CaIVC_deg = -135;              % inlet valve closing (as Simple.m)
eo.DefaultIntakePressure_bar = 1.00;
eo.DefaultIntakeTemperature_K = 293.15;
eo.DefaultRPM = 1500;

% Heat-release settings (as Simple.m)
eo.SGWindow_points = 21;          % smoothdata 'sgolay' window per cycle
eo.CombustionWindow_deg = [-20 80];
eo.SOCThreshold = 0.05;           % SOC = aROHR rises above 5% of its peak

% Exhaust analyser and molar masses [kg/mol]
eo.AnalyzerDryBasis = false;      % true if the analyser reports dry values
eo.AirO2MoleFraction = 0.2095;
eo.M_air = 28.965e-3;
eo.M_CO2 = 44.009e-3;
eo.M_NOx = 46.0055e-3;            % NOx reported as NO2 (usual convention)
eo.M_CO = 28.010e-3;
eo.M_HC = 13.88e-3;               % HC as C1 equivalent (CH1.85); check analyser

% Model validation and injection-timing prediction (faValidateModel)
eo.RunValidation = true;
eo.CalibrationRunIDs = {};            % e.g. {'B7_3bar_SOI18_R1'}; empty = per
eo.CalibrationSOI_bTDC_deg = 18;      %   group, the run nearest this SOI
eo.PredictionSOI_bTDC_deg = 14:1:21;  % timings to predict (lab range)
eo.PredictionLoadMode = 'constantLoad'; % keep gross IMEP (lab); or 'constantFuel'
eo.LoadMatchTolerance_bar = 0.3;      % same load if no TargetIMEP is given
eo.WiebeA = 6.908;                    % 99.9 % burned at the end of the Wiebe
eo.CompressionMatchOffset_deg = 5;    % match pressure 5 deg before SOI/SOC
eo.ValidationWindow_deg = [-60 120];  % pressure RMSE window
eo.MaxPeakPressure_bar = Inf;         % set if the lab gives a limit
eo.MaxPressureRise_bar_deg = Inf;     % set for a noise/pressure-rise limit

eo.RunListFile = fullfile(experimentRoot, 'Data', 'ExperimentRuns.csv');
eo.MakePlots = true;
eo.ExportResults = true;

%% Read the run list and the fuel data
runs = faReadRunList(eo.RunListFile);
expFuels = faFuels();
expFuelNames = {expFuels.Name};

%% Analyse each run
nRuns = numel(runs);
expRows = cell(nRuns, 1);
expCurves = cell(nRuns, 1);

for iRun = 1:nRuns
    thisRun = runs(iRun);

    % Fuel lookup (B7 is the 'Diesel' entry in faFuels)
    fuelName = strtrim(thisRun.Fuel);
    if strcmpi(fuelName, 'B7')
        fuelName = 'Diesel';
    end
    iFuel = find(strcmpi(expFuelNames, fuelName), 1);
    assert(~isempty(iFuel), ...
        'Run %s: fuel "%s" not in faFuels (%s).', ...
        thisRun.RunID, thisRun.Fuel, strjoin(expFuelNames, ', '));

    % Data file: relative to Data/ unless a full path is given
    thisRun.DataFilePath = thisRun.DataFile;
    if ~isfile(thisRun.DataFilePath)
        thisRun.DataFilePath = fullfile(experimentRoot, 'Data', thisRun.DataFile);
    end
    assert(isfile(thisRun.DataFilePath), 'Run %s: data file not found: %s', ...
        thisRun.RunID, thisRun.DataFile);

    fprintf('\n=== Run %d/%d: %s (%s, SOI %g deg bTDC) ===\n', ...
        iRun, nRuns, thisRun.RunID, expFuels(iFuel).Name, thisRun.SOI_bTDC_deg);

    result = faAnalyzeExperiment(thisRun, ExpCyl, SpS, expFuels(iFuel), eo);

    expRows{iRun} = result.KPI;
    expCurves{iRun} = result.Curves;
end

%% Assemble results
ExperimentResults = struct();
ExperimentResults.KPI = struct2table(vertcat(expRows{:}));
ExperimentResults.Curves = vertcat(expCurves{:});
ExperimentResults.Options = eo;
ExperimentResults.Geometry = ExpCyl;
ExperimentResults.Fuels = expFuels;

ExperimentKPI = ExperimentResults.KPI;
disp(ExperimentKPI(:, {'RunID', 'Fuel', 'SOI_bTDC_deg', 'NetIMEP_bar', ...
    'CA50_deg', 'IndicatedEfficiency_pct', 'ISFC_gkWh', 'ISCO2_gkWh'}));

%% Export (new time-stamped folder each time)
if eo.ExportResults
    outputDirectory = fullfile(experimentRoot, 'ExperimentResults', ...
        datestr(now, 'yyyymmdd_HHMMSS'));
    mkdir(outputDirectory);

    writetable(ExperimentKPI, fullfile(outputDirectory, 'Experiment_KPI.csv'));
    save(fullfile(outputDirectory, 'ExperimentResults.mat'), 'ExperimentResults');

    fprintf('\nResults saved in:\n%s\n', outputDirectory);
end

%% Plots
if eo.MakePlots
    faPlotExperiment(ExperimentResults);
end

%% Model validation and injection-timing prediction (ADDED 25/11, part 2)
% Calibrates faSimulate on one measured timing per fuel/engine/load group,
% predicts the other timings, and compares with any other measured timings.
if eo.RunValidation
    ModelValidation = faValidateModel(ExperimentResults, runs, ExpCyl, ...
        expFuels, SpS(2:5), eo);

    disp(ModelValidation.Calibration);
    disp(ModelValidation.Prediction);
    if height(ModelValidation.Validation) > 0
        disp(ModelValidation.Validation);
    end

    if eo.ExportResults
        save(fullfile(outputDirectory, 'ModelValidation.mat'), 'ModelValidation');
        writetable(ModelValidation.Calibration, ...
            fullfile(outputDirectory, 'Validation_Calibration.csv'));
        writetable(ModelValidation.Prediction, ...
            fullfile(outputDirectory, 'Validation_Prediction.csv'));
        if height(ModelValidation.Validation) > 0
            writetable(ModelValidation.Validation, ...
                fullfile(outputDirectory, 'Validation_vsMeasured.csv'));
        end
    end

    if eo.MakePlots
        faPlotValidation(ModelValidation);
    end
end
