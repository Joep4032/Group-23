%% Additional theoretical fuel analysis
% Existing project files and engine geometry are not changed.

if ~exist('Cyl', 'var') || ~isstruct(Cyl)
    error(['Define Cyl first by running the engine geometry ', ...
           'section of Simple.m.']);
end

fuelAnalysisRoot = fileparts(mfilename('fullpath'));

addpath(fullfile(fuelAnalysisRoot, 'Functions'));

FuelAnalysisOptions = faOptions();

% Optional changes:
% FuelAnalysisOptions.ComparisonMode = 'equalMass';
% FuelAnalysisOptions.MakePlots = false;
% FuelAnalysisOptions.ExportResults = false;

FuelAnalysisResults = faRunAnalysis(Cyl, FuelAnalysisOptions);
FuelAnalysisKPI = FuelAnalysisResults.KPI;

disp(FuelAnalysisKPI);

fprintf('\nTheoretical gross indicated results only.\n');
fprintf('Brake-specific quantities are unavailable: NaN.\n');
fprintf('No target IMEP: load changes with fuel and timing.\n');

if FuelAnalysisOptions.ExportResults
    fprintf('Results saved in:\n%s\n', ...
        FuelAnalysisResults.OutputDirectory);
end