cpp_path = '~/github/luxurious-loki/Verification/Output/';
matlab_path = '~/github/LoKI/Code/Output/';

step_size = 1000;
n = 1000;

disp(n);

matlab_eedf = importdata([matlab_path 'ee/' num2str(n) '/eedf.txt']);
energy = matlab_eedf.data(:,1);
matlab_eedf = matlab_eedf.data(:, 2);

cpp_eedf = importdata([cpp_path 'ee/' num2str(n) '/eedf.txt']);
cpp_eedf = cpp_eedf.data(:, 2);

rel_error = abs(matlab_eedf - cpp_eedf) ./ matlab_eedf;

figure(1);
p = plot(energy, rel_error);
set(p, 'LineWidth', 2);
set(gca, 'FontSize', 22);
xlabel('Energy (eV)');
ylabel('Relative Absolute Error');
% xlim([n n_stop] ./ 10^3);