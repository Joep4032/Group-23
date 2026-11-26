function fig = faPlotExperiment(results)
%FAPLOTEXPERIMENT Tabbed window comparing measured runs.
% ADDED 25/11 for the experimental-results analysis.
%
% Curves: one line per run. Trends: one line per group
% (fuel | target load | engine), plotted against SOI.

KPI = results.KPI;
curves = results.Curves;
nRun = numel(curves);

runLabels = cellstr(KPI.RunID);
runColors = lines(max(nRun, 1));

fig = figure( ...
    'Name', 'Experimental results', ...
    'NumberTitle', 'off', ...
    'Color', 'w', ...
    'Position', [100 60 1250 750]);

tabs = uitabgroup(fig);

%% Tab 1: mean filtered cylinder pressure

tab = uitab(tabs, 'Title', 'Pressure');
ax = axes('Parent', tab);
hold(ax, 'on');

for i = 1:nRun
    c = curves(i);
    plot(ax, c.CA, c.PressureFiltered_Pa/1e5, ...
        'Color', runColors(i,:), 'LineWidth', 1.3);
end

xlim(ax, [-60 100]);
xlabel(ax, 'Crank angle [deg aTDC]');
ylabel(ax, 'Pressure [bar absolute]');
title(ax, 'Pegged, cycle-averaged, filtered pressure');
grid(ax, 'on');
legend(ax, runLabels, 'Location', 'best', 'Interpreter', 'none');

%% Tab 2: pV diagrams (log-log)

tab = uitab(tabs, 'Title', 'pV diagrams');
ax = axes('Parent', tab);
hold(ax, 'on');

for i = 1:nRun
    c = curves(i);
    plot(ax, c.Volume_m3*1000, c.PressureFiltered_Pa/1e5, ...
        'Color', runColors(i,:), 'LineWidth', 1.3);
end

set(ax, 'XScale', 'log', 'YScale', 'log');
xlabel(ax, 'Volume [liter]');
ylabel(ax, 'Pressure [bar absolute]');
title(ax, 'Logarithmic pV diagrams');
grid(ax, 'on');
legend(ax, runLabels, 'Location', 'best', 'Interpreter', 'none');

%% Tab 3: heat release

tab = uitab(tabs, 'Title', 'Heat release');
ax1 = subplot(2,1,1, 'Parent', tab);
ax2 = subplot(2,1,2, 'Parent', tab);
hold(ax1, 'on');
hold(ax2, 'on');

for i = 1:nRun
    c = curves(i);
    plot(ax1, c.CA, c.aROHR_Jdeg, 'Color', runColors(i,:), 'LineWidth', 1.3);
    plot(ax2, c.CAwindow, c.Qnorm, 'Color', runColors(i,:), 'LineWidth', 1.3);
end

xlim(ax1, [-45 135]);
xlabel(ax1, 'Crank angle [deg aTDC]');
ylabel(ax1, 'aROHR [J/deg]');
title(ax1, 'Apparent rate of heat release (cycle-averaged, dynamic \gamma)');
grid(ax1, 'on');
legend(ax1, runLabels, 'Location', 'best', 'Interpreter', 'none');

xlabel(ax2, 'Crank angle [deg aTDC]');
ylabel(ax2, 'Normalised cumulative heat [-]');
title(ax2, 'Normalised cumulative heat release');
grid(ax2, 'on');

%% Tab 4: gamma and temperature

tab = uitab(tabs, 'Title', 'Gamma and temperature');
ax1 = subplot(2,1,1, 'Parent', tab);
ax2 = subplot(2,1,2, 'Parent', tab);
hold(ax1, 'on');
hold(ax2, 'on');

for i = 1:nRun
    c = curves(i);
    plot(ax1, c.CA, c.Gamma, 'Color', runColors(i,:), 'LineWidth', 1.3);
    plot(ax2, c.CA, c.Temperature_K, 'Color', runColors(i,:), 'LineWidth', 1.3);
end

xlim(ax1, [-180 180]);
xlim(ax2, [-180 180]);
xlabel(ax1, 'Crank angle [deg aTDC]');
ylabel(ax1, '\gamma [-]');
title(ax1, 'Dynamic heat capacity ratio (calcGamma)');
grid(ax1, 'on');
legend(ax1, runLabels, 'Location', 'best', 'Interpreter', 'none');

xlabel(ax2, 'Crank angle [deg aTDC]');
ylabel(ax2, 'Temperature [K]');
title(ax2, 'Mean in-cylinder temperature (ideal gas, from IVC)');
grid(ax2, 'on');

%% Tab 5: KPI trends against SOI, one line per group

tab = uitab(tabs, 'Title', 'SOI trends');

groupKey = KPI.Fuel + " | " + compose("%.1f bar", KPI.TargetIMEP_bar) + ...
    " | engine " + KPI.Engine;
groups = unique(groupKey, 'stable');
groupColors = lines(numel(groups));

fields = {'CA50_deg', 'IndicatedEfficiency_pct', 'ISFC_gkWh', ...
          'ISCO2_gkWh', 'ISNOx_gkWh', 'COV_IMEP_pct'};
labels = {'CA50 [deg aTDC]', 'Indicated efficiency [%]', ...
          'ISFC [g/kWh]', 'ISCO_2 [g/kWh]', 'ISNO_x [g/kWh]', ...
          'COV of net IMEP [%]'};

for iPlot = 1:numel(fields)
    ax = subplot(2,3,iPlot, 'Parent', tab);
    hold(ax, 'on');

    for g = 1:numel(groups)
        rows = sortrows(KPI(groupKey == groups(g), :), 'SOI_bTDC_deg');
        plot(ax, rows.SOI_bTDC_deg, rows.(fields{iPlot}), '-o', ...
            'Color', groupColors(g,:), 'LineWidth', 1.3);
    end

    xlabel(ax, 'SOI [deg bTDC]');
    ylabel(ax, labels{iPlot});
    grid(ax, 'on');

    if iPlot == 1
        legend(ax, cellstr(groups), 'Location', 'best', 'Interpreter', 'none');
    end
end

%% Tab 6: KPI table

tab = uitab(tabs, 'Title', 'KPI table');

tableData = table2cell(KPI);
for col = 1:width(KPI)
    if isstring(KPI{:,col})
        tableData(:,col) = cellstr(KPI{:,col});
    end
end

uitable( ...
    'Parent', tab, ...
    'Units', 'normalized', ...
    'Position', [0.01 0.06 0.98 0.93], ...
    'Data', tableData, ...
    'ColumnName', KPI.Properties.VariableNames, ...
    'RowName', []);

uicontrol( ...
    'Parent', tab, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.01 0.005 0.98 0.045], ...
    'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'left', ...
    'String', ['IS = indicated specific (pV work); BS = brake specific ', ...
    '(needs BrakePower_kW). NaN = input missing in ExperimentRuns.csv. ', ...
    'CO2 = combustion CO2 (carbon balance), not lifecycle CO2.']);
end
