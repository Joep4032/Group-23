function v = faValidateModel(expResults, runs, Cyl, fuels, species, eo)
%FAVALIDATEMODEL Calibrate the theoretical model on a measured run,
% validate it, and predict other injection timings.
% ADDED 25/11 (model validation / timing prediction).
%
% For each group of measured runs (same fuel | engine | load):
%   1. CALIBRATE faSimulate on one measured timing (the run closest to
%      eo.CalibrationSOI_bTDC_deg, or the RunIDs in eo.CalibrationRunIDs):
%        - ignition delay  = measured SOC - SOI
%        - Wiebe shape     = burn duration and exponent m fitted to the
%                            measured CA10 / CA50 / CA90
%        - compression     = initial pressure scaled so the model matches
%                            the measured pressure just before combustion
%                            (lumps heat loss / blow-by during compression)
%        - heat release    = measured apparent heat release; with fuel flow
%                            this gives etaHR = Q_apparent / (m_fuel * LHV),
%                            which lumps wall heat loss and incomplete burn
%   2. VALIDATE at that timing: measured vs model pressure, IMEP, peak
%      pressure, CA50 (shows how well the calibrated model reproduces it).
%   3. PREDICT the timings in eo.PredictionSOI_bTDC_deg with the same
%      calibration. 'constantLoad' re-solves the fuel mass so the gross
%      IMEP stays equal to the calibration point (as in the lab).
%   4. If other timings of the same group were measured, compare the
%      predictions with those measurements (true validation).
%
% Assumptions: ignition delay [ms], burn shape, etaHR and compression
% scaling stay constant over the sweep. NOx/soot are not modelled, so the
% "best" timing is judged on efficiency / fuel use and the pressure limits.
%
% Uses faSimulate, faOptions and CylinderVolume unchanged.

KPI = table2struct(expResults.KPI);
curves = expResults.Curves;
nRun = numel(KPI);

Vd = pi*Cyl.Bore^2*Cyl.Stroke/4;
fuelNames = {fuels.Name};

%% 1. Group the runs: fuel | engine | load

groupId = zeros(nRun, 1);
nGroup = 0;

for i = 1:nRun
    for g = 1:nGroup
        j = find(groupId == g, 1);
        if sameGroup(KPI(i), KPI(j), eo)
            groupId(i) = g;
            break
        end
    end
    if groupId(i) == 0
        nGroup = nGroup + 1;
        groupId(i) = nGroup;
    end
end

%% 2. Calibrate, validate and predict per group

calRows = {};
predRows = {};
valRows = {};
models = {};

for g = 1:nGroup

    members = find(groupId == g);
    calIndex = chooseCalibrationRuns(KPI, members, eo);

    for ic = calIndex(:)'

        k = KPI(ic);
        fuel = fuels(strcmp(fuelNames, char(k.Fuel)));

        m = calibrateRun(k, curves(ic), runs(ic), fuel, Cyl, species, eo);

        if isempty(m)
            continue
        end

        % Calibration-point validation
        cal = struct();
        cal.CalibrationRunID = char(k.RunID);
        cal.Fuel = char(k.Fuel);
        cal.Engine = char(k.Engine);
        cal.SOI_bTDC_deg = k.SOI_bTDC_deg;
        cal.IgnitionDelay_ms = m.Delay_s*1000;
        cal.WiebeDuration_deg = m.Options.BurnDuration_deg;
        cal.WiebeM = m.Options.WiebeM;
        cal.WiebeFitRMS_deg = m.WiebeFitRMS_deg;
        cal.HeatReleaseFraction = m.EtaHR;
        cal.HeatBasis = m.HeatBasis;
        cal.CompressionScale = m.CompressionScale;
        cal.GrossIMEP_meas_bar = k.GrossIMEP_bar;
        cal.GrossIMEP_model_bar = m.Sim.Work_J/Vd/1e5;
        cal.PeakPressure_meas_bar = k.PeakPressure_bar;
        cal.PeakPressure_model_bar = max(m.Sim.Pressure_Pa)/1e5;
        cal.CA_PeakPressure_meas_deg = k.CA_PeakPressure_deg;
        cal.CA_PeakPressure_model_deg = m.CApeakModel;
        cal.CA50_meas_deg = k.CA50_deg;
        cal.CA50_model_deg = m.Sim.CA50_degATDC;
        cal.PressureRMSE_bar = m.PressureRMSE_bar;
        cal.GrossEff_meas_pct = 100*k.GrossIMEP_bar*1e5*Vd/(m.FuelMass_kg*fuel.LHV_Jkg);
        cal.GrossEff_model_pct = 100*m.Sim.Work_J/(m.FuelMass_kg*fuel.LHV_Jkg);
        if ~m.HasFuelFlow
            cal.GrossEff_meas_pct = NaN;   % relative to apparent heat only
        end
        calRows{end+1} = cal; %#ok<AGROW>

        % Predictions over the timing sweep (+ every measured timing)
        measuredSOI = [KPI(members).SOI_bTDC_deg];
        sweep = unique([eo.PredictionSOI_bTDC_deg(:); measuredSOI(:)])';
        targetIMEP = m.Sim.Work_J/Vd/1e5;

        pred = cell(numel(sweep), 1);
        for s = 1:numel(sweep)
            pred{s} = predictTiming(m, sweep(s), targetIMEP, fuel, Cyl, ...
                species, Vd, eo);
            pred{s}.CalibrationRunID = char(k.RunID);
            pred{s}.Fuel = char(k.Fuel);
            pred{s}.Engine = char(k.Engine);
        end
        pred = vertcat(pred{:});

        % Best timing: highest efficiency that respects the limits
        ok = [pred.WithinLimits];
        eff = [pred.GrossEff_pct];
        eff(~ok) = -Inf;
        [bestEff, iBest] = max(eff);
        for s = 1:numel(pred)
            pred(s).IsBest = (s == iBest) && isfinite(bestEff);
        end

        if isfinite(bestEff)
            fprintf(['\n[%s | engine %s] calibrated on %s (SOI %g): ', ...
                'best predicted SOI = %g deg bTDC (gross eff %.1f %%, ', ...
                'peak %.1f bar)\n'], char(k.Fuel), char(k.Engine), ...
                char(k.RunID), k.SOI_bTDC_deg, pred(iBest).SOI_bTDC_deg, ...
                pred(iBest).GrossEff_pct, pred(iBest).PeakPressure_bar);
        else
            fprintf('\n[%s] no predicted timing within the pressure limits.\n', ...
                char(k.Fuel));
        end

        predRows{end+1} = pred; %#ok<AGROW>

        % Comparison with the other measured timings of this group
        for j = members(:)'
            if j == ic
                continue
            end
            kj = KPI(j);
            pj = pred([pred.SOI_bTDC_deg] == kj.SOI_bTDC_deg);
            val = struct();
            val.CalibrationRunID = char(k.RunID);
            val.MeasuredRunID = char(kj.RunID);
            val.Fuel = char(kj.Fuel);
            val.SOI_bTDC_deg = kj.SOI_bTDC_deg;
            val.CA50_meas_deg = kj.CA50_deg;
            val.CA50_pred_deg = pj.CA50_deg;
            val.PeakPressure_meas_bar = kj.PeakPressure_bar;
            val.PeakPressure_pred_bar = pj.PeakPressure_bar;
            mfj = kj.FuelFlow_gs/1000/(kj.RPM/120);
            val.FuelMass_meas_mg = mfj*1e6;
            val.FuelMass_pred_mg = pj.FuelMass_mg;
            val.GrossEff_meas_pct = 100*kj.GrossIMEP_bar*1e5*Vd/(mfj*fuel.LHV_Jkg);
            val.GrossEff_pred_pct = pj.GrossEff_pct;
            val.GrossIMEP_meas_bar = kj.GrossIMEP_bar;
            val.GrossIMEP_pred_bar = pj.GrossIMEP_bar;
            val.CA50_error_deg = val.CA50_pred_deg - val.CA50_meas_deg;
            val.PeakPressure_error_pct = 100*(val.PeakPressure_pred_bar ...
                - val.PeakPressure_meas_bar)/val.PeakPressure_meas_bar;
            val.GrossEff_error_pctpt = val.GrossEff_pred_pct - val.GrossEff_meas_pct;
            valRows{end+1} = val; %#ok<AGROW>
        end

        % Keep curves for the plots
        mm.CalibrationRunID = char(k.RunID);
        mm.Label = sprintf('%s (cal. %s, SOI %g)', char(k.Fuel), ...
            char(k.RunID), k.SOI_bTDC_deg);
        mm.Measured = curves(ic);
        mm.ModelCycle = m.Sim;
        mm.Prediction = pred;
        mm.MeasuredSOI = measuredSOI;
        mm.MeasuredRuns = members;
        models{end+1} = mm; %#ok<AGROW>
    end
end

%% 3. Assemble outputs

v.Calibration = toTable(calRows);
v.Prediction = toTable(cellfun(@(p) p(:), predRows, 'UniformOutput', false));
v.Validation = toTable(valRows);
v.Models = models;
v.MeasuredKPI = KPI;
v.Fuels = fuels;
v.Vd_m3 = Vd;
v.Options = eo;

fprintf(['\nNote: predictions keep ignition delay [ms], burn shape and ', ...
    'heat-release fraction from the calibration run.\nNOx and soot are ', ...
    'not modelled: earlier timing usually raises NOx.\n']);
end

%% ===== Local functions ====================================================

function tf = sameGroup(a, b, eo)
tf = strcmpi(char(a.Fuel), char(b.Fuel)) && ...
     strcmpi(char(a.Engine), char(b.Engine));
if ~tf
    return
end
if isfinite(a.TargetIMEP_bar) && isfinite(b.TargetIMEP_bar)
    tf = abs(a.TargetIMEP_bar - b.TargetIMEP_bar) < 1e-9;
else
    tf = abs(a.NetIMEP_bar - b.NetIMEP_bar) <= eo.LoadMatchTolerance_bar;
end
end

function idx = chooseCalibrationRuns(KPI, members, eo)
if ~isempty(eo.CalibrationRunIDs)
    ids = arrayfun(@(k) char(k.RunID), KPI(members), 'UniformOutput', false);
    idx = members(ismember(ids, eo.CalibrationRunIDs));
else
    soi = [KPI(members).SOI_bTDC_deg];
    [~, i] = min(abs(soi - eo.CalibrationSOI_bTDC_deg));
    idx = members(i);
end
end

function m = calibrateRun(k, c, run, fuel, Cyl, species, eo)
m = [];

rpm = k.RPM;
soi = k.SOI_bTDC_deg;
CaSOI = -soi;
SOC = k.SOC_deg;
target = [k.CA10_deg k.CA50_deg k.CA90_deg];

if any(~isfinite([SOC target]))
    warning('%s: SOC or CA10/50/90 missing; run not used for calibration.', ...
        char(k.RunID));
    return
end

% Ignition delay (measured)
delayDeg = max(SOC - CaSOI, 0);
delay_s = delayDeg/(6*rpm);
SOCmodel = CaSOI + delayDeg;

% Wiebe fit to CA10/50/90 (a fixed, fit m and duration; log keeps them > 0)
a = eo.WiebeA;
xs = [0.1 0.5 0.9];
zx = @(mm) (-log(1 - xs*(-expm1(-a)))/a).^(1/(mm+1));
cost = @(q) sum((SOCmodel + exp(q(2))*zx(exp(q(1))) - target).^2);
q0 = [log(2) log(max(target(3) - SOCmodel, 5)*1.3)];
q = fminsearch(cost, q0, optimset('Display', 'off', ...
    'TolX', 1e-8, 'TolFun', 1e-10, 'MaxFunEvals', 4000, 'MaxIter', 4000));
wiebeM = exp(q(1));
duration = min(exp(q(2)), 179 - SOCmodel);

% Measured apparent heat release over the combustion window
window = c.CA >= eo.CombustionWindow_deg(1) & c.CA <= eo.CombustionWindow_deg(2);
Qapp = max(cumtrapz(c.CA(window), c.aROHR_Jdeg(window)));

cyclesPerSecond = rpm/120;
hasFuelFlow = isfinite(k.FuelFlow_gs) && k.FuelFlow_gs > 0;

if hasFuelFlow
    mf = k.FuelFlow_gs/1000/cyclesPerSecond;
    etaHR = Qapp/(mf*fuel.LHV_Jkg);
    basis = 'fuel flow';
    if etaHR > 1
        warning(['%s: apparent heat release exceeds fuel energy ', ...
            '(etaHR = %.2f). Check fuel flow / pegging.'], char(k.RunID), etaHR);
    end
else
    mf = Qapp/fuel.LHV_Jkg;
    etaHR = 1;
    basis = 'apparent heat (no fuel flow)';
end

% Model options: start from faOptions and override with measured values
o = faOptions();
o.RPM = rpm;
o.ComparisonMode = 'equalMass';
o.ReferenceFuelMass_kg = mf;
o.ReferenceIgnitionDelay_s = max(delay_s, 1e-9);
o.ReferenceCetaneNumber = fuel.CetaneNumber;   % -> delay exactly as measured
o.CetaneDelayExponent = 1;
o.BurnDuration_deg = duration;
o.WiebeA = a;
o.WiebeM = wiebeM;

T_IVC = run.IntakeTemperature_K;
if ~isfinite(T_IVC)
    T_IVC = eo.DefaultIntakeTemperature_K;
end
o.T_BDC_K = T_IVC;
o.P_BDC_Pa = interp1(c.CA, c.PressureFiltered_Pa, -180);

% Heat release = etaHR x fuel energy (via an effective LHV; the fuel mass,
% and so the air/fuel ratio and CO2, stay the real ones)
fuelHR = fuel;
fuelHR.LHV_Jkg = fuel.LHV_Jkg*etaHR;

% Match the compression pressure just before combustion
thetaRef = min(SOCmodel, CaSOI) - eo.CompressionMatchOffset_deg;
pMeasRef = interp1(c.CA, c.PressureFiltered_Pa, thetaRef);
pBDC0 = o.P_BDC_Pa;

for it = 1:4
    sim = faSimulate(Cyl, fuelHR, soi, o, species);
    pModelRef = interp1(sim.Theta_degATDC, sim.Pressure_Pa, thetaRef);
    o.P_BDC_Pa = o.P_BDC_Pa*pMeasRef/pModelRef;
end

sim = faSimulate(Cyl, fuelHR, soi, o, species);

% Pressure RMSE over the validation window
inWin = c.CA >= max(eo.ValidationWindow_deg(1), -180) & ...
        c.CA <= min(eo.ValidationWindow_deg(2), 180);
pModelOnMeas = interp1(sim.Theta_degATDC, sim.Pressure_Pa, c.CA(inWin));
rmse = sqrt(mean((pModelOnMeas - c.PressureFiltered_Pa(inWin)).^2))/1e5;

[~, iMax] = max(sim.Pressure_Pa);

m.Options = o;
m.Fuel = fuel;
m.FuelHR = fuelHR;
m.FuelMass_kg = mf;
m.EtaHR = etaHR;
m.HeatBasis = basis;
m.HasFuelFlow = hasFuelFlow;
m.Delay_s = delay_s;
m.WiebeFitRMS_deg = sqrt(cost(q)/3);
m.CompressionScale = o.P_BDC_Pa/pBDC0;
m.Sim = sim;
m.CApeakModel = sim.Theta_degATDC(iMax);
m.PressureRMSE_bar = rmse;
end

function p = predictTiming(m, soi, targetIMEP, fuel, Cyl, species, Vd, eo)
o = m.Options;

if strcmpi(eo.PredictionLoadMode, 'constantLoad')
    % Re-solve the fuel mass so the gross IMEP equals the calibration point
    imepAt = @(mf) simulateIMEP(mf, soi, m.FuelHR, o, Cyl, species, Vd) - targetIMEP;

    % Oxygen limit (same idea as faRunLoadStudy)
    V_BDC = CylinderVolume(Cyl.TDCangle - 180, Cyl);
    nO2 = o.P_BDC_Pa*V_BDC/(8.314*o.T_BDC_K)*o.AirO2MoleFraction;
    atomicMass = [12.011 1.008 15.999]*1e-3;
    fuelMolarMass = sum(fuel.CHO.*atomicMass);
    oxygenDemand = fuel.CHO(1) + fuel.CHO(2)/4 - fuel.CHO(3)/2;
    maxMass = 0.98*nO2/oxygenDemand*fuelMolarMass;

    lo = 0.2*m.FuelMass_kg;
    hi = min(2*m.FuelMass_kg, maxMass);
    while imepAt(hi) < 0 && hi < maxMass
        hi = min(1.5*hi, maxMass);
    end
    assert(imepAt(lo) < 0 && imepAt(hi) >= 0, ...
        'Cannot hold %.2f bar gross IMEP at SOI %g within the oxygen limit.', ...
        targetIMEP, soi);
    mf = fzero(imepAt, [lo hi], optimset('Display', 'off', 'TolX', 1e-13));
else
    mf = m.FuelMass_kg;   % constantFuel
end

o.ReferenceFuelMass_kg = mf;
sim = faSimulate(Cyl, m.FuelHR, soi, o, species);

dpdtheta = gradient(sim.Pressure_Pa, sim.Theta_degATDC)/1e5;

p.CalibrationRunID = '';
p.Fuel = '';
p.Engine = '';
p.SOI_bTDC_deg = soi;
p.LoadMode = eo.PredictionLoadMode;
p.FuelMass_mg = mf*1e6;
p.GrossIMEP_bar = sim.Work_J/Vd/1e5;
p.GrossEff_pct = 100*sim.Work_J/(mf*fuel.LHV_Jkg);
p.GrossISFC_gkWh = mf/sim.Work_J*3.6e9;
p.GrossISCO2_gkWh = sim.CO2Mass_kg/sim.Work_J*3.6e9;
p.CA50_deg = sim.CA50_degATDC;
p.PeakPressure_bar = max(sim.Pressure_Pa)/1e5;
p.MaxPressureRise_bar_deg = max(dpdtheta);
p.Lambda = sim.Lambda;
p.WithinLimits = p.PeakPressure_bar <= eo.MaxPeakPressure_bar && ...
                 p.MaxPressureRise_bar_deg <= eo.MaxPressureRise_bar_deg;
p.IsBest = false;
end

function imep = simulateIMEP(mf, soi, fuelHR, o, Cyl, species, Vd)
o.ReferenceFuelMass_kg = mf;
sim = faSimulate(Cyl, fuelHR, soi, o, species);
imep = sim.Work_J/Vd/1e5;
end

function T = toTable(rows)
% Cell array of structs (or struct arrays) -> table, one row per struct.
if isempty(rows)
    T = table();
    return
end
T = struct2table(vertcat(rows{:}), 'AsArray', true);
end
