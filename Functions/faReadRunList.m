function runs = faReadRunList(file)
%FAREADRUNLIST Read the experiment run list (Data/ExperimentRuns.csv).
% ADDED 25/11 for the experimental-results analysis.
%
% One row per measured run (repeats get their own row). Empty numeric
% cells become NaN; the related KPIs are then reported as NaN.
% Rows with Include = 0 are skipped (e.g. runs with bad data).

assert(isfile(file), 'Run list not found: %s', file);

opts = detectImportOptions(file, 'FileType', 'text', 'Delimiter', ',');

textColumns = {'RunID', 'DataFile', 'Fuel', 'Engine', 'Notes'};
textColumns = textColumns(ismember(textColumns, opts.VariableNames));
numericColumns = setdiff(opts.VariableNames, textColumns, 'stable');

opts = setvartype(opts, textColumns, 'char');
opts = setvartype(opts, numericColumns, 'double');

T = readtable(file, opts);

required = {'RunID', 'DataFile', 'Fuel', 'SOI_bTDC_deg'};
missing = required(~ismember(required, T.Properties.VariableNames));
assert(isempty(missing), 'Run list is missing column(s): %s', ...
    strjoin(missing, ', '));

% Optional columns that the analysis reads
optionalNumeric = {'Include', 'RPM', 'TargetIMEP_bar', ...
    'IntakePressure_bar', 'IntakeTemperature_K', 'FuelFlow_gs', ...
    'AirFlow_gs', 'CO2_volpct', 'CO_ppm', 'HC_ppm', 'NOx_ppm', ...
    'BrakePower_kW'};

for k = 1:numel(optionalNumeric)
    if ~ismember(optionalNumeric{k}, T.Properties.VariableNames)
        T.(optionalNumeric{k}) = NaN(height(T), 1);
    end
end

for k = {'Engine', 'Notes'}
    if ~ismember(k{1}, T.Properties.VariableNames)
        T.(k{1}) = repmat({''}, height(T), 1);
    end
end

% Skip excluded rows and rows without a data file
keep = ~cellfun(@isempty, T.DataFile) & ~(T.Include == 0);
T = T(keep, :);

assert(height(T) > 0, 'No runs to analyse in %s.', file);

runs = table2struct(T);
end
