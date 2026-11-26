function c = faSimulate(Cyl, fuel, soiBTDC, o, species)
%FASIMULATE Closed compression/expansion cycle with prescribed combustion.
%
% Integrates from compression BDC (-180 deg) to expansion BDC (+180 deg).
% Constant-volume closure contributes zero work.
%
% Work and efficiency are idealized GROSS INDICATED quantities.
% No target IMEP or brake/friction model is used.

global Runiv

%% Fuel mass

if strcmp(o.ComparisonMode, 'equalEnergy')
    energyReference = o.ReferenceFuelMass_kg*o.ReferenceLHV_Jkg;
    mf = energyReference/fuel.LHV_Jkg;
else
    mf = o.ReferenceFuelMass_kg;
end

fuelEnergy = mf*fuel.LHV_Jkg;

%% Geometry: reuse original CylinderVolume function

theta = (-180:o.Step_deg:180)';

% theta is relative to compression TDC.
% CylinderVolume expects the angle corresponding to Cyl.TDCangle.
V = CylinderVolume(theta+Cyl.TDCangle, Cyl);

assert(all(isfinite(V) & V > 0), 'Invalid cylinder volumes.');

%% Stoichiometry and initial air charge

M = [species.Mass]; % kg/mol: O2, N2, CO2, H2O

% Atomic masses consistent with the supplied NASA database.
MO = M(1)/2;
MH = (M(4)-MO)/2;
MC = M(3)-2*MO;

carbon = fuel.CHO(1);
hydrogen = fuel.CHO(2);
oxygen = fuel.CHO(3);

fuelMolarMass = carbon*MC + hydrogen*MH + oxygen*MO;
oxygenDemand = carbon + hydrogen/4 - oxygen/2;

assert(fuelMolarMass > 0 && oxygenDemand > 0 && ...
    fuel.LHV_Jkg > 0 && fuel.CetaneNumber > 0, ...
    'Invalid fuel properties.');

nAir = o.P_BDC_Pa*V(1)/(Runiv*o.T_BDC_K);

initialMass = nAir * ...
    [o.AirO2MoleFraction, 1-o.AirO2MoleFraction, 0, 0] .* M;

nFuel = mf/fuelMolarMass;

massChange = nFuel * ...
    [-oxygenDemand, 0, carbon, hydrogen/2] .* M;

lambda = (initialMass(1)/M(1))/(nFuel*oxygenDemand);

assert(lambda >= 1, ...
    'Complete combustion requires lambda >= 1. Reduce the fuel input.');

assert(abs(sum(massChange)-mf) < 1e-10*mf, ...
    'Stoichiometric mass balance failed.');

%% Assumed ignition delay and prescribed combustion

delay_s = o.ReferenceIgnitionDelay_s * ...
    (o.ReferenceCetaneNumber/fuel.CetaneNumber)^o.CetaneDelayExponent;

SOC = -soiBTDC + 6*o.RPM*delay_s;

assert(SOC > -180 && SOC+o.BurnDuration_deg < 180, ...
    'The complete burn must lie inside the simulated cycle.');

[xb, dxb] = faWiebe(theta, SOC, o.BurnDuration_deg, ...
    o.WiebeA, o.WiebeM);

%% Thermodynamic model inputs

model.Options = o;
model.Species = species;
model.InitialGasMass_kg = initialMass;
model.GasMassChange_kg = massChange;
model.FuelMass_kg = mf;

model.ReferenceEnthalpy_Jkg = zeros(1,4);

for k = 1:4
    model.ReferenceEnthalpy_Jkg(k) = ...
        HNasa(o.ReferenceTemperature_K, species(k));
end

%% Initial state

n = numel(theta);

T = zeros(n,1);
p = zeros(n,1);
U = zeros(n,1);
gamma = zeros(n,1);
work = zeros(n,1);

T(1) = o.T_BDC_K;

state = faMixture(T(1), xb(1), model);

U(1) = state.U_J;
p(1) = state.nR_JK*T(1)/V(1);
gamma(1) = state.Gamma;

%% Energy-conserving crank-angle integration

for k = 2:n

    dV = V(k)-V(k-1);
    dQ = fuelEnergy*(xb(k)-xb(k-1));

    trialT = T(k-1);
    converged = false;

    for iteration = 1:o.MaximumIterations

        state = faMixture(trialT, xb(k), model);
        trialP = state.nR_JK*trialT/V(k);

        residual = state.U_J-U(k-1) ...
            + 0.5*(p(k-1)+trialP)*dV-dQ;

        if abs(residual) <= o.EnergyTolerance_J
            converged = true;
            break
        end

        derivative = state.Cv_JK ...
            + 0.5*state.nR_JK/V(k)*dV;

        assert(isfinite(derivative) && derivative > 0, ...
            'Invalid energy-equation derivative; reduce Step_deg.');

        trialT = trialT-residual/derivative;

        assert(isfinite(trialT) && ...
            trialT >= o.TemperatureBounds_K(1) && ...
            trialT <= o.TemperatureBounds_K(2), ...
            ['Temperature left the configured interval. ', ...
             'Check inputs or reduce Step_deg.']);
    end

    assert(converged, ...
        'Energy solver did not converge at CA = %.3f degrees.', theta(k));

    T(k) = trialT;
    p(k) = trialP;
    U(k) = state.U_J;
    gamma(k) = state.Gamma;

    work(k) = work(k-1) + 0.5*(p(k-1)+p(k))*dV;
end

%% Combustion phasing and diagnostics

z50 = (-log(1-0.5*(-expm1(-o.WiebeA)))/o.WiebeA) ...
    ^(1/(o.WiebeM+1));

CA50 = SOC + o.BurnDuration_deg*z50;

energyResidual = U(end)-U(1)+work(end) ...
    - fuelEnergy*(xb(end)-xb(1));

%% Results

c.FuelName = fuel.Name;
c.SOI_bTDC_deg = soiBTDC;

c.Theta_degATDC = theta;
c.Volume_m3 = V;
c.Pressure_Pa = p;
c.Temperature_K = T;
c.Gamma = gamma;

c.BurnFraction = xb;
c.HeatRelease_Jdeg = fuelEnergy*dxb;
c.CumulativeHeat_J = fuelEnergy*xb;

c.InternalEnergy_J = U;
c.CumulativeWork_J = work;
c.Work_J = work(end);

c.FuelMass_kg = mf;
c.FuelEnergy_J = fuelEnergy;
c.CO2Mass_kg = massChange(3);
c.AirMass_kg = sum(initialMass);
c.Lambda = lambda;

c.SOC_degATDC = SOC;
c.CA50_degATDC = CA50;
c.IgnitionDelay_ms = delay_s*1000;

c.EnergyBalanceError_pct = ...
    100*abs(energyResidual)/fuelEnergy;
end