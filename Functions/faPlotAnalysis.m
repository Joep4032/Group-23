function fig = faPlotAnalysis(r)
%FAPLOTANALYSIS Display the additional analysis in one tabbed window.

o = r.Options;
fuels = r.Fuels;
cycles = r.Cycles;
KPI = r.KPI;

nFuel = numel(fuels);
nTiming = numel(o.SOI_bTDC_deg);

fuelColors = lines(nFuel);
timingColors = lines(nTiming);

fuelLabels = {fuels.Name};

timingLabels = arrayfun( ...
    @(x) sprintf('%g deg BTDC', x), ...
    o.SOI_bTDC_deg, 'UniformOutput', false);

referenceIndex = find( ...
    abs(o.SOI_bTDC_deg-o.ReferenceSOI_bTDC_deg) < 1e-9, 1);

maxPressure = 0;
maxVolume = 0;
minVolume = Inf;

for i = 1:nFuel
    for j = 1:nTiming
        c = cycles{i,j};

        maxPressure = max(maxPressure, max(c.Pressure_Pa)/1e5);
        maxVolume = max(maxVolume, max(c.Volume_m3)*1000);
        minVolume = min(minVolume, min(c.Volume_m3)*1000);
    end
end

pressureLimit = ceil(maxPressure/10)*10;
volumeLimit = ceil(maxVolume*10)/10;

fig = figure( ...
    'Name', 'Theoretical fuel and injection-timing analysis', ...
    'NumberTitle', 'off', ...
    'Color', 'w', ...
    'Position', [80 80 1250 750]);

tabs = uitabgroup(fig);

%% Tab 1: fuel comparison pV diagrams

tab = uitab(tabs, 'Title', 'Fuel pV comparison');

ax1 = subplot(1,2,1, 'Parent', tab);
ax2 = subplot(1,2,2, 'Parent', tab);

hold(ax1, 'on');
hold(ax2, 'on');

for i = 1:nFuel
    c = cycles{i,referenceIndex};

    % Append constant-volume idealized closure.
    % This is not an intake/exhaust or pumping-loop model.
    volume = [c.Volume_m3; c.Volume_m3(1)]*1000;
    pressure = [c.Pressure_Pa; c.Pressure_Pa(1)]/1e5;

    plot(ax1, volume, pressure, ...
        'Color', fuelColors(i,:), 'LineWidth', 1.5);

    plot(ax2, volume, pressure, ...
        'Color', fuelColors(i,:), 'LineWidth', 1.5);
end

set(ax2, 'XScale', 'log', 'YScale', 'log');

xlabel(ax1, 'Volume [liter]');
ylabel(ax1, 'Pressure [bar absolute]');
title(ax1, 'Linear pV Diagram');

xlabel(ax2, 'Volume [liter]');
ylabel(ax2, 'Pressure [bar absolute]');
title(ax2, 'Logarithmic pV diagram');

xlim(ax1, [0 volumeLimit]);
ylim(ax1, [0 pressureLimit]);

xlim(ax2, [0.9*minVolume volumeLimit]);
ylim(ax2, [0.5 pressureLimit]);

grid(ax1, 'on');
grid(ax2, 'on');

legend(ax1, fuelLabels, 'Location', 'best');
legend(ax2, fuelLabels, 'Location', 'best');

% %% Tab 2: pressure and prescribed heat release
% 
% tab = uitab(tabs, 'Title', 'Pressure and heat release');
% 
% ax1 = subplot(2,1,1, 'Parent', tab);
% ax2 = subplot(2,1,2, 'Parent', tab);
% 
% hold(ax1, 'on');
% hold(ax2, 'on');
% 
% for i = 1:nFuel
%     c = cycles{i,referenceIndex};
% 
%     plot(ax1, c.Theta_degATDC, c.Pressure_Pa/1e5, ...
%         'Color', fuelColors(i,:), 'LineWidth', 1.5);
% 
%     plot(ax2, c.Theta_degATDC, c.HeatRelease_Jdeg, ...
%         'Color', fuelColors(i,:), 'LineWidth', 1.5);
% end
% 
% xlabel(ax1, 'Crank angle relative to compression TDC [deg]');
% ylabel(ax1, 'Pressure [bar absolute]');
% title(ax1, sprintf('SOI %g deg BTDC', o.ReferenceSOI_bTDC_deg));
% 
% xlabel(ax2, 'Crank angle relative to compression TDC [deg]');
% ylabel(ax2, 'Prescribed heat release [J/deg]');
% title(ax2, 'Wiebe heat release; not measured aROHR');
% 
% xlim(ax1, [-40 100]);
% xlim(ax2, [-40 100]);
% 
% grid(ax1, 'on');
% grid(ax2, 'on');
% 
% legend(ax1, fuelLabels, 'Location', 'best');
% legend(ax2, fuelLabels, 'Location', 'best');


%% Tab 3: temperature and dynamic gamma

tab = uitab(tabs, 'Title', 'Temperature and gamma');

ax1 = subplot(2,1,1, 'Parent', tab);
ax2 = subplot(2,1,2, 'Parent', tab);

hold(ax1, 'on');
hold(ax2, 'on');

for i = 1:nFuel
    c = cycles{i,referenceIndex};

    plot(ax1, c.Theta_degATDC, c.Temperature_K, ...
        'Color', fuelColors(i,:), 'LineWidth', 1.5);

    plot(ax2, c.Theta_degATDC, c.Gamma, ...
        'Color', fuelColors(i,:), 'LineWidth', 1.5);
end

xlabel(ax1, 'Crank angle relative to compression TDC [deg]');
ylabel(ax1, 'Temperature [K]');
title(ax1, 'Temperature from the energy equation');

xlabel(ax2, 'Crank angle relative to compression TDC [deg]');
ylabel(ax2, 'Gas-phase \gamma [-]');
title(ax2, 'NASA temperature and composition-dependent properties');

grid(ax1, 'on');
grid(ax2, 'on');

legend(ax1, fuelLabels, 'Location', 'best');
legend(ax2, fuelLabels, 'Location', 'best');

%% Tab 4: timing pV diagrams

tab = uitab(tabs, 'Title', 'SOI pV sweep');

for i = 1:nFuel
    % CHANGED 25/11: grid grows with the number of fuels (was 2x2, which
    % fails once the two B7-OME3 blends are added: 6 fuels -> 2x3).
    ax = subplot(2, ceil(nFuel/2), i, 'Parent', tab);
    hold(ax, 'on');

    for j = 1:nTiming
        c = cycles{i,j};

        volume = [c.Volume_m3; c.Volume_m3(1)]*1000;
        pressure = [c.Pressure_Pa; c.Pressure_Pa(1)]/1e5;

        plot(ax, volume, pressure, ...
            'Color', timingColors(j,:), 'LineWidth', 1.4);
    end

    xlabel(ax, 'Volume [liter]');
    ylabel(ax, 'Pressure [bar absolute]');
    title(ax, fuels(i).Name);

    xlim(ax, [0 volumeLimit]);
    ylim(ax, [0 pressureLimit]);

    grid(ax, 'on');
    legend(ax, timingLabels, 'Location', 'best');
end

%% Tab 5: timing effects on work and peak pressure

tab = uitab(tabs, 'Title', 'Timing trends');

ax1 = subplot(1,2,1, 'Parent', tab);
ax2 = subplot(1,2,2, 'Parent', tab);

hold(ax1, 'on');
hold(ax2, 'on');

for i = 1:nFuel

    index = (i-1)*nTiming + (1:nTiming);

    plot(ax1, -o.SOI_bTDC_deg, KPI.GrossWork_JCycle(index), ...
        '-o', 'Color', fuelColors(i,:), 'LineWidth', 1.5);

    plot(ax2, -o.SOI_bTDC_deg, KPI.PeakPressure_bar(index), ...
        '-o', 'Color', fuelColors(i,:), 'LineWidth', 1.5);
end

xlabel(ax1, 'SOI [deg relative to TDC; negative = BTDC]');
ylabel(ax1, 'Gross indicated work [J/cycle]');
title(ax1, 'Effect of injection timing on gross indicated work');

xlabel(ax2, 'SOI [deg relative to TDC; negative = BTDC]');
ylabel(ax2, 'Peak pressure [bar absolute]');
title(ax2, 'Effect of injection timing on peak cylinder pressure');

grid(ax1, 'on');
grid(ax2, 'on');

legend(ax1, fuelLabels, 'Location', 'best');
legend(ax2, fuelLabels, 'Location', 'best');

%% Tab 6: KPI table

tab = uitab(tabs, 'Title', 'KPI results');

tableData = table2cell(KPI);
tableData(:,1) = cellstr(KPI.Fuel);

uitable( ...
    'Parent', tab, ...
    'Units', 'normalized', ...
    'Position', [0.01 0.05 0.98 0.94], ...
    'Data', tableData, ...
    'ColumnName', KPI.Properties.VariableNames, ...
    'RowName', [], ...
    'ColumnWidth', 'auto');

uicontrol( ...
    'Parent', tab, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.01 0.005 0.98 0.035], ...
    'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'left', ...
    'String', ['NaN = unavailable brake quantity. ', ...
    'CO2 is complete-combustion exhaust CO2, not lifecycle emissions.']);
end