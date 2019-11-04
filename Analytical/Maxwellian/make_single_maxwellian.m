% --- DATA FOR THE TESTS ---
k0 = 1e-12;                   % elastic momentum transfer constant rate coefficient
Tg = 300;                     % gas temperature to be considered in the calculations
Muam = 40;                    % mass of the heavy species to be considered in the calculations (in UAM)
decades = 15;                 % decades of decay considered for the eedf
reducedField = (0:1:50);      % range of reduced electric fields to be considered in the calculations

% --- CONSTANTS USED IN THE TESTS ---
kb = 1.38064852e-23;
me = 9.10938356e-31;
e = 1.6021766208e-19;
uam = 1.660539040e-27;

% --- THEORETICAL CALCULATIONS ---
M = Muam*uam;
Te = kb*Tg/e+M/(3*me^2*e)*(e*reducedField*1e-21/k0).^2;
TeSimulation = zeros(size(Te));

cellNumber = 2000;

for i = 1:length(reducedField)
        
        % --- EVALUATE ENERGY GRID OF THE SIMULATION ---
        maxEnergy = str2double(sprintf('%.2e', (decades/log10(exp(1)))*Te(i)));
        energy = linspace(0, maxEnergy, cellNumber+1);
        energy(cellNumber+2) = 2*energy(end)-energy(end-1);
        
        % --- READ SIMULATION RESULTS ---
        outputFolderStr = sprintf('../../build/Output/Maxwellian%scellNumber_%d%sreducedField_%f', filesep, cellNumber, ...
            filesep, reducedField(i));
        
        fid = fopen([outputFolderStr filesep 'swarm_parameters.txt'], 'r');
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        info = strsplit(fgetl(fid));
        TeSimulation(i) = str2double(info{5});
        fclose(fid);
        
        fid = fopen([outputFolderStr filesep 'eedf.txt'], 'r');
        fgetl(fid);
        data = fscanf(fid,'%e',[3 inf]);
        theoreticalEEDF = exp(-data(1,:)/Te(i))/sum(sqrt(data(1,:)).*exp(-data(1,:)/Te(i)))/(data(1,2)-data(1,1));
        relErrorEEDF(i) = sum(abs(theoreticalEEDF-data(2,:))./theoreticalEEDF)/cellNumber;
        fclose(fid);
end

clf;
set(gcf, 'DefaultAxesFontSize', 18);
fig = figure(1);

black = [0 0 0];
set(fig,'defaultAxesColorOrder',[black; black]);

yyaxis left;
ax_one = plot(reducedField, relErrorEEDF .* 1e4);
xlabel('Reduced Electric Field (Td)');
ylabel('EEDF Average Relative Error (\times10^{-4})');

yyaxis right;
ax_two = plot(reducedField, abs(TeSimulation - Te)./Te .* 1e4, '--');
ylabel('Electron Temperature Relative Error (\times10^{-4})');

% set(ax_one,'color',	[0, 0.4470, 0.7410]);
set(ax_one, 'LineWidth', 3);
set(ax_two,'color', 'r');%[0.8500, 0.3250, 0.0980]);
set(ax_two, 'LineWidth', 3);

ylim([1.783 1.81]);