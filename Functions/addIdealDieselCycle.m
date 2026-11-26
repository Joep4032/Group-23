function addIdealDieselCycle(ax, Cyl)
%ADDIDEALDIESELCYCLE Overlay an ideal Diesel cycle on existing pV axes.
%
% Existing axes must use:
%   Volume: dm^3 or litres
%   Pressure: bar absolute
%
% Benchmark assumptions:
%   Initial pressure = 1 bar absolute
%   Constant gamma = 1.35
%   Cut-off ratio V3/V2 = 2
%
% This is an idealized reference, not a fit to the measured cycle.

%% Assumptions
p1 = 1e5;                 % Pa absolute
gammaIdeal = 1.35;
cutoffRatio = 2.0;

r = Cyl.CompressionRatio;

assert(cutoffRatio > 1 && cutoffRatio < r, ...
    'Cut-off ratio must be between 1 and the compression ratio.');

%% Reuse existing engine geometry
Vd = pi*Cyl.Bore^2*Cyl.Stroke/4;
Vc = Vd/(r-1);

V1 = Vc+Vd;               % Start of compression
V2 = Vc;                  % End of compression
V3 = cutoffRatio*V2;      % End of heat addition
V4 = V1;                  % End of expansion

%% State pressures
p2 = p1*r^gammaIdeal;
p3 = p2;
p4 = p3*(V3/V4)^gammaIdeal;

%% Compression and expansion trajectories
Vcompression = linspace(V1, V2, 500);
pCompression = p1*(V1./Vcompression).^gammaIdeal;

Vexpansion = linspace(V3, V4, 500);
pExpansion = p3*(V3./Vexpansion).^gammaIdeal;

%% Overlay without deleting existing plots
wasHeld = ishold(ax);
hold(ax, 'on');

hCompression = plot(ax, ...
    Vcompression*1000, pCompression/1e5, ...
    'Color', [0 0.75 0], 'LineWidth', 1.8);

hExpansion = plot(ax, ...
    Vexpansion*1000, pExpansion/1e5, ...
    'Color', [0 0.75 0], 'LineWidth', 1.8);

hAddition = plot(ax, ...
    [V2 V3]*1000, [p2 p3]/1e5, ...
    'r-', 'LineWidth', 2);

hRejection = plot(ax, ...
    [V4 V1]*1000, [p4 p1]/1e5, ...
    'c-', 'LineWidth', 2);

% Optional idealized intake/exhaust reference:
% Both strokes have the same pressure, so net pumping work is zero.
hGasExchange = plot(ax, ...
    [V2 V1]*1000, [p1 p1]/1e5, ...
    'm-', 'LineWidth', 1.5);

%% Legend
% Identify the original measured line, excluding the new overlays.
allLines = findobj(ax, 'Type', 'line');

idealLines = [
    hCompression
    hExpansion
    hAddition
    hRejection
    hGasExchange
];

measuredLines = allLines(~ismember(allLines, idealLines));

if ~isempty(measuredLines)
    legend(ax, ...
        [measuredLines(1), hCompression, ...
         hAddition, hRejection, hGasExchange], ...
        {'Measured cycle', ...
         'Ideal compression/expansion', ...
         'Ideal heat addition', ...
         'Ideal heat rejection', ...
         'Ideal intake/exhaust reference'}, ...
        'Location', 'best', 'FontSize', 8);
end

%% Expand limits so the ideal pressure peak is visible
limits = ylim(ax);

upperPressure = 10*ceil(max(limits(2), p2/1e5)/10);

if strcmp(ax.YScale, 'log')
    ylim(ax, [0.5 upperPressure]);
else
    ylim(ax, [0 upperPressure]);
end

xlabel(ax, 'Volume [liter]');
ylabel(ax, 'Pressure [bar absolute]');
grid(ax, 'on');

if ~wasHeld
    hold(ax, 'off');
end
end