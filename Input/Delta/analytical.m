moM = 2.5e-5;
sigma0 = 1e-19;
Q = 1e-19;
T = 300;
EoN = 10;
uMax = 30;
n = 2999;
U0 = 10;

e = 1.6021766208e-19;

step = uMax / n;
W = e / sqrt(moM * 6) * EoN / Q;

%Vector::LinSpaced(cellNumber, .5, cellNumber - .5) * step


u = linspace(.5, n - .5, n) * step;

g = 1/(moM*4)*sigma0/Q*U0/W^2*exp(-U0/(2*W^2)) .* (expint(u.^2 ./ (2 * W^2)) - expint(U0^2 / (2 * W^2)));

%to be set by normalization
fref = 0;
eref = 0;

fH0 = fref * exp(-U0^2/(2*W^2) * ((u./U0).^2 - (eref/U0)^2));

%steps to solving this equation :
% 1. choose an initial value for eref
% 2. normalize to find fref
% 3. solve fref = f0 to find new value for eref
% 4. iterate until abs(f0 - f0_old)/f0 < a threshold.
f0 = fH0 .* (1 + g);