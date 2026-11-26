function s = faMixture(T, xb, model)
%FAMIXTURE Thermodynamic state of reacting air/products.
%
% Gas species order: O2, N2, CO2, H2O.
% Unburned fuel is a pressureless thermal reservoir.
%
% uGas = h(T)-h(Tref)-R*T.
% The chemical-energy source is the specified fuel LHV.
%
% This simplified model excludes vaporization, spray mixing,
% dissociation, wall heat transfer and gas exchange.

global Runiv

bounds = model.Options.TemperatureBounds_K;

assert(isfinite(T) && T >= bounds(1) && T <= bounds(2), ...
    'Temperature is outside the configured NASA evaluation interval.');

mass = model.InitialGasMass_kg + xb*model.GasMassChange_kg;

assert(all(mass >= -1e-14), ...
    'Negative species mass: check oxygen availability.');

cp = zeros(1,4);
cv = zeros(1,4);
u = zeros(1,4);

Rspecies = Runiv ./ [model.Species.Mass];

for k = 1:4
    cp(k) = CpNasa(T, model.Species(k));
    cv(k) = CvNasa(T, model.Species(k));

    u(k) = HNasa(T, model.Species(k)) ...
        - model.ReferenceEnthalpy_Jkg(k) ...
        - Rspecies(k)*T;
end

unburnedMass = model.FuelMass_kg*(1-xb);
fuelCv = model.Options.UnburnedFuelHeatCapacity_JkgK;

fuelU = fuelCv*(T-model.Options.ReferenceTemperature_K);

s.U_J = sum(mass.*u) + unburnedMass*fuelU;

s.Cv_JK = sum(mass.*cv) + unburnedMass*fuelCv;

s.nR_JK = sum(mass.*Rspecies);

% Gas-phase heat-capacity ratio, reported as a diagnostic.
s.Gamma = sum(mass.*cp)/sum(mass.*cv);

assert(s.Cv_JK > 0 && s.nR_JK > 0, ...
    'Invalid mixture heat capacity or gas constant.');
end