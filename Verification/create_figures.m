cpp_path = '~/github/luxurious-loki/Verification/Output/';
matlab_path = '~/github/LoKI/Code/Output/';

step_size = 1000;
n_start = 1000;
n_stop = 8000;

N = n_start:step_size:n_stop;

ave_rel_error = zeros(size(N));

disp(size(ave_rel_error));

for n = N
    disp(n);

    matlab_eedf = importdata([matlab_path 'growth/' num2str(n) '/eedf.txt']);
    matlab_eedf = matlab_eedf.data(:, 2);

    cpp_eedf = importdata([cpp_path 'growth/' num2str(n) '/eedf.txt']);
    cpp_eedf = cpp_eedf.data(:, 2);

    rel_error = abs(matlab_eedf - cpp_eedf) ./ matlab_eedf;

    ave_rel_error((n-n_start) / step_size + 1) = sum(rel_error) / n;
end

set(gcf, 'DefaultAxesFontSize', 18);
figure(1);
p = plot(N ./ 10^3, ave_rel_error);
xlabel('Grid size (\times10^{3})');
ylabel('Average Relative Absolute Error');
xlim([n_start n_stop] ./ 10^3);