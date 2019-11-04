n = 10;
step = 0.01;
offset = 9.95;

x = offset + (0:n)*step;
y = zeros(1, n+1);

y(n/2+1) = 1e-18;

set(gcf, 'DefaultAxesFontSize', 18);
figure(1);
clf;
plot(x, y);
xlabel('Energy (eV)');
ylabel('Distribution (eV^{-1})');

% u_max = 30
% n = 3000
% delta_u = .01

x_cells = offset + ((1:n)-.5)*step;
y_cells = zeros(1, n);

y_cells(n/2) = 5e-19;
y_cells(n/2+1) = 5e-19;

figure(2);
clf;
plot(x_cells, y_cells);
xlabel('Energy (eV)');
ylabel('Distribution (eV^{-1})');
ylim([0, 1e-18]);