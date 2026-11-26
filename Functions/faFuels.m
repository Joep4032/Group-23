function fuels = faFuels()
%FAFUELS Representative properties for theoretical comparisons.
%
% Diesel, FT diesel and HVO are mixtures. Their CHO entries are
% surrogate formulas used for oxygen demand and carbon balance.
%
% Pure OME3:
% CH3-O-(CH2-O)3-CH3 = C5H12O4.
%
% Density is included for reference but is not used to set fuel mass.
% Update these inputs when verified fuel specifications are available.

names = {'Diesel', 'FT diesel', 'OME3', 'HVO'};

LHV = [43.0 44.1 19.6 44.0] * 1e6;  % J/kg
CN = [52.2 75 73 75];
density = [835 780 1040 780];          % kg/m^3

% Columns: carbon, hydrogen, oxygen atoms per molecule.
CHO = [
    14 28 0
    16 34 0
     5 12 4
    16 34 0
];

fuels = struct([]);

for k = 1:numel(names)
    fuels(k).Name = names{k};
    fuels(k).LHV_Jkg = LHV(k);
    fuels(k).CetaneNumber = CN(k);
    fuels(k).Density_kgm3 = density(k);
    fuels(k).CHO = CHO(k,:);
end

% FT diesel is a theoretical comparison fuel, not a supplied test fuel.

%% ===== ADDED 25/11: B7/OME3 blends tested in the lab =====================
% The lab blends are made with commercial B7 (Q&A 20/11). The 'Diesel'
% entry above already uses B7 handbook values (43 MJ/kg, CN 52.2), so it is
% used as the base fuel. Blend ratios are in VOLUME percent of OME3.
%
% Mixing rules (assumptions, ideal mixing):
%   mass fraction  w_OME = x*rho_OME / (x*rho_OME + (1-x)*rho_B7)
%   LHV            mass-weighted   (exact for ideal mixing)
%   density        volume-weighted (ideal mixing)
%   cetane number  volume-weighted (linear blending; approximation)
%   CHO            mole-weighted pseudo-molecule, so the carbon mass
%                  fraction and oxygen demand per kg match the blend.

blendVolumeFractions = [0.125 0.25];
blendNames = {'B7-OME3 12.5%', 'B7-OME3 25%'};

iBase = find(strcmp({fuels.Name}, 'Diesel'), 1);
iOME = find(strcmp({fuels.Name}, 'OME3'), 1);

for b = 1:numel(blendVolumeFractions)
    fuels(end+1) = faBlend(fuels(iBase), fuels(iOME), ...
        blendVolumeFractions(b), blendNames{b}); %#ok<AGROW>
end
%% ===== END ADDED ===========================================================
end

%% ===== ADDED 25/11: local blending function ================================
function blend = faBlend(base, additive, x, name)
% x = volume fraction of the additive (e.g. 0.125 for 12.5 vol% OME3).

rhoBase = base.Density_kgm3;
rhoAdd = additive.Density_kgm3;

% Mass fraction of the additive
wAdd = x*rhoAdd/(x*rhoAdd + (1-x)*rhoBase);
wBase = 1-wAdd;

% Atomic masses [g/mol]; only ratios matter for the pseudo-molecule
atomicMass = [12.011 1.008 15.999];
molarMassBase = sum(base.CHO.*atomicMass);
molarMassAdd = sum(additive.CHO.*atomicMass);

% Moles of each component per kg of blend -> mole fractions
nBase = wBase/molarMassBase;
nAdd = wAdd/molarMassAdd;
yAdd = nAdd/(nBase + nAdd);

blend.Name = name;
blend.LHV_Jkg = wBase*base.LHV_Jkg + wAdd*additive.LHV_Jkg;
blend.CetaneNumber = (1-x)*base.CetaneNumber + x*additive.CetaneNumber;
blend.Density_kgm3 = (1-x)*rhoBase + x*rhoAdd;
blend.CHO = (1-yAdd)*base.CHO + yAdd*additive.CHO;
end
%% ===== END ADDED ===========================================================