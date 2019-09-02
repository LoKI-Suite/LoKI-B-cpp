% initialize parameters according to section 3.5
moM = 2.5e-5;
sigma0 = 1e-19;
Q = 1e-19;
T = 300;
EoN = 10;
U0 = 10;

e = 1.6021766208e-19;

% definitions related to grid
uMax = 30;
n = 3000;

step = uMax / n;
W = e / sqrt(moM * 6) * EoN / Q;

% define the energy grid
u = linspace(.5, n - .5, n) * step;

% construct g following (3.18) with fixed exponential integral following
% matlab doc
g = 1/(moM*4) * sigma0/Q * U0^2/W^2 * exp(-U0^2/(2*W^2)) * (Ei(u.^2 ./ (2 * W^2)) - Ei(U0^2 / (2 * W^2)));

% set epsilon values according to section 3.4.1
e0 = U0;
eref = U0;

% construct fH0 for eref = U0 and e0 = U0 following (3.11)
fH0 = exp(-e0^2/(2*W^2) * ((u./e0).^2 - (eref/e0)^2));

% construct full expression
f0 = fH0 .* (1 + g);

% normalize the expression over the given domain following (3.4)
f0 = f0 ./ (dot(f0, sqrt(u)) * step);

% check that expression is indeed normalized
disp(dot(f0, sqrt(u)));
disp(dot(f0, sqrt(u)) * step);

% log-plot the result
semilogy(u, f0);