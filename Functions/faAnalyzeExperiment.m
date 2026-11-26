function r = faAnalyzeExperiment(run, Cyl, SpS, fuel, eo)
%FAANALYZEEXPERIMENT Analyse one measured engine run.
% ADDED 25/11. Uses the same steps as the measured-data part of Simple.m
% (pegging/averaging/filtering, CylinderVolume, calcGamma, aROHR,
% CA10/50/90), but for any data file and with that run's own SOI,
% intake conditions, fuel and measured flows.
%
% Inputs:
%   run   One row of the run list (faReadRunList), with DataFilePath added.
%   Cyl   Engine geometry (same struct as Simple.m).
%   SpS   NASA species {Diesel, O2, N2, CO2, H2O} (as loaded in Simple.m).
%   fuel  Fuel entry from faFuels (LHV, CHO, ...), incl. the B7-OME3 blends.
%   eo    Options set in RunExperimentAnalysis.m.
%
% Output r.KPI (one row) and r.Curves (mean traces for plotting).
%
% Indicated KPIs come from the pV work (net = full 720 deg, incl. pumping).
% Brake KPIs need the measured brake/generator power (BrakePower_kW).
% Emission mass flows from analyser concentrations need the air flow.

%% 1. Read the pressure file and split it into cycles

data = readmatrix(run.DataFilePath, 'FileType', 'text');

CAall = data(:, eo.CrankAngleColumn);
pAll = data(:, eo.PressureColumn)*eo.PressureUnit_Pa;

valid = isfinite(CAall) & isfinite(pAll);   % drops header/empty rows
CAall = CAall(valid);
pAll = pAll(valid);

d = diff(CAall);
step = median(d(d > 0));
pointsPerCycle = round(720/step);

assert(abs(pointsPerCycle*step-720) < 1e-6, ...
    '%s: crank-angle step %.4g deg does not divide 720.', run.RunID, step);

% Start at the first sample of a full cycle (e.g. -360 deg)
first = find(abs(CAall-eo.CycleStart_deg) < step/2, 1);

assert(~isempty(first), ...
    ['%s: no sample at %g deg. Check the crank-angle convention ', ...
     '(eo.CycleStart_deg).'], run.RunID, eo.CycleStart_deg);

nCycles = floor((numel(CAall)-first+1)/pointsPerCycle);

assert(nCycles >= 1, '%s: less than one full cycle in the file.', ...
    run.RunID);

used = first:(first+nCycles*pointsPerCycle-1);
droppedSamples = numel(CAall)-numel(used);

if droppedSamples > 0
    warning('%s: %d samples outside full cycles were ignored.', ...
        run.RunID, droppedSamples);
end

Ca = reshape(CAall(used), pointsPerCycle, nCycles);
p = reshape(pAll(used), pointsPerCycle, nCycles);

%% 2. Pegging, averaging and filtering (existing function)

pIntake = valueOrDefault(run.IntakePressure_bar, ...
    eo.DefaultIntakePressure_bar, 'IntakePressure_bar', run.RunID)*1e5;

processed = processMeasuredPressure(Ca, p, pIntake);

CA = processed.CA;
pPegged = processed.PeggedPressure;
pFiltered = processed.FilteredPressure;

V = CylinderVolume(CA, Cyl);
Vd = pi*Cyl.Bore^2*Cyl.Stroke/4;

%% 3. Work and IMEP per cycle

Wnet = trapz(V, pPegged);                          % full 720 deg
IMEPnet = Wnet/Vd/1e5;

grossRange = CA >= -180 & CA <= 180;               % compression + expansion
Wgross = trapz(V(grossRange), pPegged(grossRange,:));
IMEPgross = Wgross/Vd/1e5;

%% 4. Heat release with dynamic gamma (same method as Simple.m)

T_IVC = valueOrDefault(run.IntakeTemperature_K, ...
    eo.DefaultIntakeTemperature_K, 'IntakeTemperature_K', run.RunID);

CaSOI = -run.SOI_bTDC_deg;                         % negative = before TDC

dVdCA = gradient(V, CA);

totalROHR = zeros(pointsPerCycle, 1);
totalGamma = zeros(pointsPerCycle, 1);
totalT = zeros(pointsPerCycle, 1);

for k = 1:nCycles
    pCycle = smoothdata(pPegged(:,k), 'sgolay', eo.SGWindow_points);
    dpdCA = gradient(pCycle, CA);

    [gamma, T] = calcGamma(CA, pCycle, V, eo.CaIVC_deg, CaSOI, T_IVC, SpS);

    aROHR = (gamma./(gamma-1)).*pCycle.*dVdCA + ...
            (1./(gamma-1)).*V.*dpdCA;

    totalROHR = totalROHR + aROHR;
    totalGamma = totalGamma + gamma;
    totalT = totalT + T;
end

aROHR = totalROHR/nCycles;
gammaMean = totalGamma/nCycles;
TMean = totalT/nCycles;

% CA10 / CA50 / CA90 in the combustion window (as in Simple.m)
window = CA >= eo.CombustionWindow_deg(1) & CA <= eo.CombustionWindow_deg(2);
CAwin = CA(window);
Qcum = cumtrapz(CAwin, aROHR(window));

if max(Qcum) > 0
    Qnorm = Qcum./max(Qcum);
    [Qunique, iu] = unique(Qnorm);
    CA10 = interp1(Qunique, CAwin(iu), 0.10);
    CA50 = interp1(Qunique, CAwin(iu), 0.50);
    CA90 = interp1(Qunique, CAwin(iu), 0.90);
else
    Qnorm = NaN(size(Qcum));
    [CA10, CA50, CA90] = deal(NaN);
end

% Start of combustion (SOC): last crank angle between SOI and the aROHR
% peak where the aROHR is still below eo.SOCThreshold x peak. A threshold
% is used because small modelling errors can keep the aROHR slightly
% positive before the real combustion starts.
searchRange = find(CA >= CaSOI & CA <= eo.CombustionWindow_deg(2));
[peakROHR, iPeakLocal] = max(aROHR(searchRange));
iPeak = searchRange(iPeakLocal);
beforePeak = searchRange(searchRange <= iPeak);
iBelow = beforePeak(aROHR(beforePeak) < eo.SOCThreshold*peakROHR);

if isempty(iBelow) || peakROHR <= 0
    SOC = NaN;
else
    SOC = CA(min(iBelow(end)+1, iPeak));
end

%% 5. Pressure indicators (filtered mean cycle)

[pMax, iMax] = max(pFiltered);
dpFiltered = gradient(pFiltered, CA)/1e5;          % bar/deg
maxRise = max(dpFiltered(window));

%% 6. Fuel, power and efficiency

rpm = valueOrDefault(run.RPM, eo.DefaultRPM, 'RPM', run.RunID);
cyclesPerSecond = rpm/120;                         % 4-stroke

atomicMass = [12.011 1.008 15.999]*1e-3;           % kg/mol
fuelMolarMass = sum(fuel.CHO.*atomicMass);         % kg/mol (pseudo-molecule)
carbonFraction = fuel.CHO(1)*atomicMass(1)/fuelMolarMass;
oxygenDemand = fuel.CHO(1) + fuel.CHO(2)/4 - fuel.CHO(3)/2;  % mol O2/mol fuel

Wmean = mean(Wnet);
Pind = Wmean*cyclesPerSecond;                      % W

mdotFuel = run.FuelFlow_gs/1000;                   % kg/s (NaN if missing)
fuelPower = mdotFuel*fuel.LHV_Jkg;                 % W

etaInd = Pind/fuelPower;
ISFC = mdotFuel/Pind*3.6e9;                        % g/kWh

% Combustion CO2 from the carbon balance (complete combustion)
CO2perKgFuel = carbonFraction*eo.M_CO2/atomicMass(1);
mdotCO2cb = mdotFuel*CO2perKgFuel;                 % kg/s
ISCO2 = mdotCO2cb/Pind*3.6e9;
combustionCO2_gMJ = CO2perKgFuel/fuel.LHV_Jkg*1e9;

%% 7. Air flow, lambda and analyser emissions (optional inputs)

mdotAir = run.AirFlow_gs/1000;                     % kg/s
AFRstoich = oxygenDemand/eo.AirO2MoleFraction*eo.M_air/fuelMolarMass;
lambda = (mdotAir/mdotFuel)/AFRstoich;

% Exhaust molar flow (wet): air + net mole increase from combustion
nFuel = mdotFuel/fuelMolarMass;
nExhaust = mdotAir/eo.M_air + nFuel*(fuel.CHO(2)/4 + fuel.CHO(3)/2);

% Analysers often measure on a dry basis; convert to wet if requested
if eo.AnalyzerDryBasis
    xH2O = nFuel*fuel.CHO(2)/2/nExhaust;
    wetFactor = 1-xH2O;
else
    wetFactor = 1;
end

mdotCO2meas = run.CO2_volpct/100*wetFactor*nExhaust*eo.M_CO2;
mdotNOx = run.NOx_ppm*1e-6*wetFactor*nExhaust*eo.M_NOx;
mdotCO = run.CO_ppm*1e-6*wetFactor*nExhaust*eo.M_CO;
mdotHC = run.HC_ppm*1e-6*wetFactor*nExhaust*eo.M_HC;

%% 8. Brake values (only if the brake/generator power was measured)

Pbrake = run.BrakePower_kW*1000;                   % W

%% 9. Collect the KPI row

kpi = struct();
kpi.RunID = string(run.RunID);
kpi.Fuel = string(fuel.Name);
kpi.Engine = string(run.Engine);
kpi.RPM = rpm;
kpi.TargetIMEP_bar = run.TargetIMEP_bar;
kpi.SOI_bTDC_deg = run.SOI_bTDC_deg;
kpi.Cycles = nCycles;

kpi.NetIMEP_bar = mean(IMEPnet);
kpi.GrossIMEP_bar = mean(IMEPgross);
kpi.PumpingIMEP_bar = mean(IMEPnet-IMEPgross);
kpi.COV_IMEP_pct = 100*std(IMEPnet)/abs(mean(IMEPnet));

kpi.PeakPressure_bar = pMax/1e5;
kpi.CA_PeakPressure_deg = CA(iMax);
kpi.MaxPressureRise_bar_deg = maxRise;

kpi.SOC_deg = SOC;
kpi.IgnitionDelay_deg = SOC-CaSOI;
kpi.IgnitionDelay_ms = (SOC-CaSOI)/(6*rpm)*1000;
kpi.CA10_deg = CA10;
kpi.CA50_deg = CA50;
kpi.CA90_deg = CA90;
kpi.BurnDuration10to90_deg = CA90-CA10;

kpi.IndicatedPower_kW = Pind/1000;
kpi.FuelFlow_gs = run.FuelFlow_gs;
kpi.FuelPower_kW = fuelPower/1000;
kpi.IndicatedEfficiency_pct = 100*etaInd;
kpi.ISFC_gkWh = ISFC;
kpi.CombustionCO2_gMJ = combustionCO2_gMJ;
kpi.ISCO2_gkWh = ISCO2;

kpi.Lambda = lambda;
kpi.CO2Measured_vs_CarbonBalance = mdotCO2meas/mdotCO2cb;
kpi.ISNOx_gkWh = mdotNOx/Pind*3.6e9;
kpi.ISCO_gkWh = mdotCO/Pind*3.6e9;
kpi.ISHC_gkWh = mdotHC/Pind*3.6e9;

kpi.BrakePower_kW = run.BrakePower_kW;
kpi.MechanicalEfficiency_pct = 100*Pbrake/Pind;
kpi.BrakeEfficiency_pct = 100*Pbrake/fuelPower;
kpi.BSFC_gkWh = mdotFuel/Pbrake*3.6e9;
kpi.BSCO2_gkWh = mdotCO2cb/Pbrake*3.6e9;
kpi.BSNOx_gkWh = mdotNOx/Pbrake*3.6e9;
kpi.BSCO_gkWh = mdotCO/Pbrake*3.6e9;
kpi.BSHC_gkWh = mdotHC/Pbrake*3.6e9;

kpi.PeggingOffset_bar = mean(processed.PeggingOffsets_Pa)/1e5;
kpi.DataFile = string(run.DataFile);
kpi.Notes = string(run.Notes);

r.KPI = kpi;

%% 10. Mean curves for plotting

r.Curves.RunID = run.RunID;
r.Curves.CA = CA;
r.Curves.Volume_m3 = V;
r.Curves.PressureFiltered_Pa = pFiltered;
r.Curves.aROHR_Jdeg = aROHR;
r.Curves.CAwindow = CAwin;
r.Curves.Qnorm = Qnorm;
r.Curves.Gamma = gammaMean;
r.Curves.Temperature_K = TMean;
r.Curves.IMEPnet_bar = IMEPnet;
end

%% Local helper
function value = valueOrDefault(value, default, name, runID)
if ~isfinite(value)
    warning('%s: %s missing, using default %g.', runID, name, default);
    value = default;
end
end
