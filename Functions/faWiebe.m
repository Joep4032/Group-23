function [xb, dxb] = faWiebe(theta, SOC, duration, a, m)
%FAWIEBE Normalized burned fraction and heat-release shape.
% theta, SOC and duration are in crank-angle degrees.
% dxb has units of 1/degree.

z = (theta-SOC)/duration;

xb = zeros(size(theta));
dxb = zeros(size(theta));

inside = z > 0 & z < 1;
zz = z(inside);

normalizer = -expm1(-a);

xb(inside) = -expm1(-a*zz.^(m+1))/normalizer;

dxb(inside) = a*(m+1)/duration .* zz.^m .* ...
    exp(-a*zz.^(m+1))/normalizer;

xb(z >= 1) = 1;
end