function [gamma, T_cycle] = calcGamma(CA_cycle, p_cycle, V_cycle, CaIVC, CaSOI, T_IVC, SpS)
% Computes the dynamic heat capacity ratio gamma(theta) over a single engine cycle.


    % 1. Gas composition Y(theta) for [O2, N2, CO2, H2O] (mass fractions)
	% These are assumed values
    Y_air     = [0.233, 0.767, 0.000, 0.000];
    Y_exhaust = [0.100, 0.750, 0.080, 0.070];

    % Linear transition between CaSOI and 40 deg CA after TDC
    progress = (CA_cycle - CaSOI) ./ (40 - CaSOI);
    progress = max(0, min(1, progress));

    Y_cycle = (1 - progress) .* Y_air + progress .* Y_exhaust;

    % 2. Mixture molar mass M_mix(theta) from SpS(2:5) = {O2, N2, CO2, H2O}
    M_species = [SpS(2).Mass, SpS(3).Mass, SpS(4).Mass, SpS(5).Mass];
    M_mix = 1 ./ sum(Y_cycle ./ M_species, 2);

    % 3. In-cylinder temperature T(theta) via ideal gas law from IVC
    [~, idx_IVC] = min(abs(CA_cycle - CaIVC));
    p_IVC = p_cycle(idx_IVC);
    V_IVC = V_cycle(idx_IVC);
    M_IVC = M_mix(idx_IVC);

    T_cycle = T_IVC .* (p_cycle .* V_cycle ./ (p_IVC * V_IVC)) .* (M_mix ./ M_IVC);
	T_cycle = max(250, min(3500, T_cycle)); % Clamp within NASA polynomial bounds else NASA function returns NaN is T is outside range 200-3500

    % 4. Mixture specific heats cp and cv using NASA polynomials
    cp_mix = Y_cycle(:,1) .* CpNasa(T_cycle, SpS(2)) + ...
             Y_cycle(:,2) .* CpNasa(T_cycle, SpS(3)) + ...
             Y_cycle(:,3) .* CpNasa(T_cycle, SpS(4)) + ...
             Y_cycle(:,4) .* CpNasa(T_cycle, SpS(5));

    cv_mix = Y_cycle(:,1) .* CvNasa(T_cycle, SpS(2)) + ...
             Y_cycle(:,2) .* CvNasa(T_cycle, SpS(3)) + ...
             Y_cycle(:,3) .* CvNasa(T_cycle, SpS(4)) + ...
             Y_cycle(:,4) .* CvNasa(T_cycle, SpS(5));

    % 5. Dynamic heat capacity ratio
    gamma = cp_mix ./ cv_mix;
end