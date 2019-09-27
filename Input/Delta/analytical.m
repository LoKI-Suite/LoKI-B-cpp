e = 1.6021766208e-19;

eV = e;
Td = 1e-21;

% initialize parameters according to section 3.5
moM = 2.5e-5;
sigma0 = 1e-20;
Q = 1e-19;
T = 300;
EoN = 10 * Td;
U0 = 10.00 * eV;

% definitions related to grid
uMax = 30 * eV;
n = 6000;

step = uMax / n;
W = e / sqrt(moM * 6) * EoN / Q;

% define the energy grid
u = linspace(.5, n - .5, n) * step;

% construct g following (3.18) with fixed exponential integral following
% matlab doc (+ Ei terms have been switched)
nU0 = floor(U0 / step);

g = zeros(1, n);
g(1, 1:nU0-1) = 1/(moM*4) * sigma0/Q * U0^2/W^2 * exp(-U0^2/(2 * W^2)) * (Ei(U0^2 / (2 * W^2)) - Ei(u(1, 1:nU0-1).^2 / (2 * W^2)));

% set epsilon values according to section 3.4.1
e0 = U0;
eref = U0;

% construct fH0 for eref = U0 and e0 = U0 following (3.11)
fH0 = exp(-e0^2/(2*W^2) * ((u./e0).^2 - (eref/e0)^2));

% construct full expression
f0 = fH0 .* (1 + g);

% normalize the expression over the given domain following (3.4)
f0 = f0 / (dot(f0, sqrt(u)) * step);

% check that expression is indeed normalized
% disp(dot(f0, sqrt(u)));
% disp(dot(f0, sqrt(u)) * step);

% plot g
% figure;
% plot(u/eV, g);

eedf = importdata('eedf_matlab_2_6000.txt');
eedf = eedf.data;

% log-plot the result
figure;
subplot(1,2,1); semilogy(u / eV, f0 * eV^(3/2));
hold on;
semilogy(u / eV, eedf(:,2));

subplot(1,2,2);
semilogy(u/eV, abs(f0 * eV^(3/2) - eedf(:,2)') ./ (f0 * eV^(3/2)));