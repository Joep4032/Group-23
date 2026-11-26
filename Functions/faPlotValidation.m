function fig = faPlotValidation(v)
%FAPLOTVALIDATION Plots for faValidateModel.
% ADDED 25/11 (model validation / timing prediction).
%   Tab 1: calibration run, measured vs calibrated model
%   Tab 2: predicted KPIs vs SOI, measured points overlaid
%   Tab 3: tables (calibration, prediction, validation)

models = v.Models;
nModel = numel(models);

fig = figure( ...
    'Name', 'Model validation and timing prediction', ...
    'NumberTitle', 'off', ...
    'Color', 'w', ...
    'Position', [120 40 1250 780]);

tabs = uitabgroup(fig);

if nModel == 0
    tab = uitab(tabs, 'Title', 'No calibration');
    uicontrol('Parent', tab, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.05 0.45 0.9 0.1], 'BackgroundColor', 'w', ...
        'String', 'No run could be used for calibration (see warnings).');
    return
end

colors = lines(nModel);

%% Tab 1: calibration, measured vs model

tab = uitab(tabs, 'Title', 'Calibration');

for i = 1:nModel
    mm = models{i};
    c = mm.Measured;
    s = mm.ModelCycle;

    ax1 = subplot(nModel, 2, 2*i-1, 'Parent', tab);
    hold(ax1, 'on');
    plot(ax1, c.CA, c.PressureFiltered_Pa/1e5, 'k', 'LineWidth', 1.5);
    plot(ax1, s.Theta_degATDC, s.Pressure_Pa/1e5, '--', ...
        'Color', colors(i,:), 'LineWidth', 1.5);
    xlim(ax1, [-60 120]);
    xlabel(ax1, 'Crank angle [deg aTDC]');
    ylabel(ax1, 'Pressure [bar]');
    title(ax1, mm.Label, 'Interpreter', 'none');
    legend(ax1, {'Measured', 'Calibrated model'}, 'Location', 'best');
    grid(ax1, 'on');

    ax2 = subplot(nModel, 2, 2*i, 'Parent', tab);
    hold(ax2, 'on');
    plot(ax2, c.CA, c.aROHR_Jdeg, 'k', 'LineWidth', 1.5);
    plot(ax2, s.Theta_degATDC, s.HeatRelease_Jdeg, '--', ...
        'Color', colors(i,:), 'LineWidth', 1.5);
    xlim(ax2, [-30 100]);
    xlabel(ax2, 'Crank angle [deg aTDC]');
    ylabel(ax2, 'Heat release rate [J/deg]');
    title(ax2, 'Measured aROHR vs fitted Wiebe');
    legend(ax2, {'Measured aROHR', 'Model (Wiebe)'}, 'Location', 'best');
    grid(ax2, 'on');
end

%% Tab 2: timing prediction with measured points

tab = uitab(tabs, 'Title', 'Timing prediction');

predFields = {'GrossEff_pct', 'GrossISFC_gkWh', 'CA50_deg', ...
              'PeakPressure_bar', 'MaxPressureRise_bar_deg', 'FuelMass_mg'};
labels = {'Gross indicated efficiency [%]', 'Gross ISFC [g/kWh]', ...
          'CA50 [deg aTDC]', 'Peak pressure [bar]', ...
          'Max pressure rise [bar/deg]', 'Fuel mass [mg/cycle]'};

% Measured values on the same basis (NaN where not available)
K = v.MeasuredKPI;
measuredValue = @(k, field) measuredFor(k, field, v);

legendItems = {};
legendHandles = [];

for iPlot = 1:numel(predFields)
    ax = subplot(2, 3, iPlot, 'Parent', tab);
    hold(ax, 'on');

    for i = 1:nModel
        mm = models{i};
        P = mm.Prediction;
        soi = [P.SOI_bTDC_deg];
        [soi, order] = sort(soi);
        y = [P.(predFields{iPlot})];
        y = y(order);

        h = plot(ax, soi, y, '-', 'Color', colors(i,:), 'LineWidth', 1.5);

        % Best predicted timing
        best = find([P(order).IsBest], 1);
        if ~isempty(best)
            plot(ax, soi(best), y(best), 'p', 'MarkerSize', 12, ...
                'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
        end

        % Measured points of this group (calibration run = square)
        for j = mm.MeasuredRuns(:)'
            yj = measuredValue(K(j), predFields{iPlot});
            if isfinite(yj)
                if strcmp(char(K(j).RunID), mm.CalibrationRunID)
                    marker = 's';
                else
                    marker = 'o';
                end
                plot(ax, K(j).SOI_bTDC_deg, yj, marker, 'MarkerSize', 7, ...
                    'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i,:));
            end
        end

        if iPlot == 1
            legendHandles(end+1) = h; %#ok<AGROW>
            legendItems{end+1} = mm.Label; %#ok<AGROW>
        end
    end

    xlabel(ax, 'SOI [deg bTDC]');
    ylabel(ax, labels{iPlot});
    grid(ax, 'on');
end

ax = subplot(2, 3, 1, 'Parent', tab);
legend(ax, legendHandles, legendItems, 'Location', 'best', 'Interpreter', 'none');

uicontrol('Parent', tab, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.01 0.0 0.98 0.03], 'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'left', ...
    'String', ['Line = calibrated model prediction; square = calibration ', ...
    'run; circles = other measured timings; star = best predicted timing ', ...
    '(max efficiency within pressure limits). NOx not modelled.']);

%% Tab 3: tables

tab = uitab(tabs, 'Title', 'Tables');
addTable(tab, v.Calibration, [0.01 0.70 0.98 0.27], 'Calibration and validation at the calibration timing');
addTable(tab, v.Prediction, [0.01 0.36 0.98 0.30], 'Predicted timings');
addTable(tab, v.Validation, [0.01 0.02 0.98 0.30], 'Prediction vs other measured timings');
end

%% Local functions

function y = measuredFor(k, field, v)
Vd = v.Vd_m3;
switch field
    case 'GrossEff_pct'
        mf = k.FuelFlow_gs/1000/(k.RPM/120);
        y = 100*k.GrossIMEP_bar*1e5*Vd/(mf*lhvOf(k, v));
    case 'GrossISFC_gkWh'
        mf = k.FuelFlow_gs/1000/(k.RPM/120);
        y = mf/(k.GrossIMEP_bar*1e5*Vd)*3.6e9;
    case 'CA50_deg'
        y = k.CA50_deg;
    case 'PeakPressure_bar'
        y = k.PeakPressure_bar;
    case 'MaxPressureRise_bar_deg'
        y = k.MaxPressureRise_bar_deg;
    case 'FuelMass_mg'
        y = k.FuelFlow_gs/(k.RPM/120)*1000;
    otherwise
        y = NaN;
end
end

function L = lhvOf(k, v)
idx = find(strcmp({v.Fuels.Name}, char(k.Fuel)), 1);
L = NaN;
if ~isempty(idx)
    L = v.Fuels(idx).LHV_Jkg;
end
end

function addTable(parent, T, position, heading)
uicontrol('Parent', parent, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [position(1) position(2)+position(4) position(3) 0.025], ...
    'BackgroundColor', 'w', 'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold', 'String', heading);

if isempty(T) || width(T) == 0
    data = {'(none yet: measure more timings of the same fuel and load)'};
    names = {''};
else
    data = table2cell(T);
    names = T.Properties.VariableNames;
end

uitable('Parent', parent, 'Units', 'normalized', 'Position', position, ...
    'Data', data, 'ColumnName', names, 'RowName', []);
end
