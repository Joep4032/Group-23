function o = faOptions()
%FAOPTIONS Configuration for theoretical fuel and timing analysis.
% Units are SI unless explicitly stated.

o.RPM = 1500;

% Positive values mean degrees BEFORE compression TDC.
o.SOI_bTDC_deg = [14 16 18 21];
o.ReferenceSOI_bTDC_deg = 18;

% 'equalEnergy' or 'equalMass'
o.ComparisonMode = 'equalEnergy';

% Previous reference input:
% 24 mg/cycle * 43 MJ/kg = 1032 J/cycle.
o.ReferenceFuelMass_kg = 24e-6;
o.ReferenceLHV_Jkg = 43e6;

% Ideal closed-cycle initial state at compression BDC.
% These are assumptions, not measured intake conditions.
o.P_BDC_Pa = 1e5;
o.T_BDC_K = 293.15;
o.AirO2MoleFraction = 0.21;

% Assumed fuel-dependent ignition-delay scaling:
% delay = referenceDelay * (referenceCN / fuelCN)^exponent.
o.ReferenceIgnitionDelay_s = 0.8e-3;
o.ReferenceCetaneNumber = 52.2;
o.CetaneDelayExponent = 1;

% Prescribed combustion parameters.
o.BurnDuration_deg = 35;
o.WiebeA = 6.908;
o.WiebeM = 2;

% Thermal reference for sensible enthalpy calculations.
o.ReferenceTemperature_K = 298.15;

% Unburned fuel is represented by a pressureless thermal reservoir.
% This generic value is an assumption, not a measured fuel property.
% Injection hydraulics and evaporation are not modelled.
o.UnburnedFuelHeatCapacity_JkgK = 2000;

% Numerical settings.
o.Step_deg = 0.2;
o.EnergyTolerance_J = 1e-8;
o.MaximumIterations = 30;
o.TemperatureBounds_K = [200 3500];

o.MakePlots = true;
o.ExportResults = true;
end