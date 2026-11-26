function result = processMeasuredPressure(Ca, p, referencePressure)
%PROCESSMEASUREDPRESSURE Peg, average and filter measured pressure.
%
% Inputs:
%   Ca                Crank-angle matrix [deg], one cycle per column.
%   p                 Pressure matrix [Pa], one cycle per column.
%   referencePressure Assumed/measured intake pressure [Pa absolute].
%
% Processing:
%   1. Peg each cycle using the mean pressure around intake BDC.
%   2. Average all pegged cycles.
%   3. Apply a 21-point cubic Savitzky-Golay filter.
%
% The default pegging window is -182 to -178 degrees.
% A 1 bar reference is an assumption unless independently measured.

validateattributes(referencePressure, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});

assert(isequal(size(Ca), size(p)), ...
    'Pressure and crank-angle matrices must have the same size.');

assert(all(isfinite(Ca(:))) && all(isfinite(p(:))), ...
    'Input data contain NaN or Inf values.');

CA = Ca(:,1);

assert(all(diff(CA) > 0), ...
    'Crank angles must increase within each cycle.');

gridDifference = bsxfun(@minus, Ca, CA);

assert(max(abs(gridDifference(:))) < 1e-8, ...
    'Cycle crank-angle grids differ. Resample them before averaging.');

step = median(diff(CA));

assert(max(abs(diff(CA)-step)) < 1e-8, ...
    'Savitzky-Golay filtering requires a uniform crank-angle grid.');

%% Peg each individual cycle
peggingWindow = abs(CA+180) <= 2;

assert(any(peggingWindow), ...
    'No data found in the pegging window around -180 degrees.');

offsets = referencePressure-mean(p(peggingWindow,:), 1);

pPegged = bsxfun(@plus, p, offsets);

%% Ensemble-average pressure
pMean = mean(pPegged, 2);

%% Cubic Savitzky-Golay filtering: 21 points
windowLength = 21;
polynomialOrder = 3;
halfWindow = (windowLength-1)/2;

nPoints = numel(CA);

assert(nPoints >= windowLength, ...
    'There are too few samples for the 21-point filter.');

pFiltered = zeros(size(pMean));

% Fit a local polynomial and evaluate it at the current sample.
% At the endpoints, use asymmetric windows rather than zero padding.
for i = 1:nPoints

    firstIndex = max(1, min(i-halfWindow, ...
        nPoints-windowLength+1));

    indices = (firstIndex:firstIndex+windowLength-1)';

    localPosition = indices-i;

    A = ones(windowLength, polynomialOrder+1);

    for order = 1:polynomialOrder
        A(:,order+1) = localPosition.^order;
    end

    coefficients = A\pMean(indices);

    pFiltered(i) = coefficients(1);
end

assert(all(pFiltered > 0), ...
    ['Processed absolute pressure is non-positive. ', ...
     'Check pressure units and the pegging reference.']);

%% Return results
result.CA = CA;
result.RawPressure = p;
result.PeggedPressure = pPegged;
result.MeanPressure = pMean;
result.FilteredPressure = pFiltered;
result.PeggingOffsets_Pa = offsets;

result.ReferencePressure_Pa = referencePressure;
result.NumberOfCycles = size(p,2);
result.FilterWindow_points = windowLength;
result.FilterOrder = polynomialOrder;

fprintf('\nMeasured pressure processing\n');
fprintf('Cycles averaged: %d\n', result.NumberOfCycles);
fprintf('Pegging reference: %.3f bar absolute\n', ...
    referencePressure/1e5);
fprintf('Mean pegging offset: %.3f bar\n', mean(offsets)/1e5);
fprintf('Filter: cubic Savitzky-Golay, %d points\n', windowLength);
fprintf('Filter first-to-last sample span: %.2f deg CA\n', ...
    (windowLength-1)*step);

end