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

result = zeros(100, 2);

for cellNumber = 100:100:10000
    %     mkdir(sprintf('Maxwellian/cellNumber_%i', cellNumber));
    
    disp(cellNumber);
    
    for i = 1:length(reducedField)
        
        % --- EVALUATE ENERGY GRID OF THE SIMULATION ---
        maxEnergy = str2double(sprintf('%.2e', (decades/log10(exp(1)))*Te(i)));
        energy = linspace(0, maxEnergy, cellNumber+1);
        energy(cellNumber+2) = 2*energy(end)-energy(end-1);
        
        % --- READ SIMULATION RESULTS ---
        outputFolderStr = sprintf('../../build/Output/Maxwellian%scellNumber_%d%sreducedField_%f', filesep, cellNumber, ...
            filesep, reducedField(i));
        
        fid = fopen([outputFolderStr '/swarm_parameters.txt'], 'r');
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
    
%     summaryFolder = sprintf('Maxwellian_summary/cellNumber_%i', cellNumber);
%     
%     mkdir(summaryFolder);
%     
%     fid = fopen([summaryFolder filesep 'summary.txt'], 'w');
%     fprintf(fid, '################################\n');
%     fprintf(fid, '# DATA CONSIDERED FOR THE TESTS#\n');
%     fprintf(fid, '################################\n\n');
%     fprintf(fid, '# elastic momentum transfer constant rate coefficient = %e (m3s-1)\n', k0);
%     fprintf(fid, '#                                     gas temperature = %f (K)\n', Tg);
%     fprintf(fid, '#                           mass of the heavy species = %f (uam)\n', Muam);
%     fprintf(fid, '#               decades of decay ensured for the eedf = %d \n', decades);
%     fprintf(fid, '#                 number of cells for the energy grid = %d \n\n', cellNumber);
%     fprintf(fid, 'E/N(Td)              Te_th(eV)            Te_sim(eV)           RelError(Te)         RelError(EEDF)\n');
%     values(5:5:5*length(Te)) = relErrorEEDF;
%     values(4:5:5*length(Te)) = abs(Te-TeSimulation)./Te;
%     values(3:5:5*length(Te)) = TeSimulation;
%     values(2:5:5*length(Te)) = Te;
%     values(1:5:5*length(Te)) = reducedField;
%     fprintf(fid, '%#.14e %#.14e %#.14e %#.14e %#.14e \n', values);
%     fclose(fid);
%     hold on;
%     plot(reducedField, (Te-TeSimulation)./Te, reducedField, relErrorEEDF);
    
    result(cellNumber/100, 1) = sum(relErrorEEDF) / length(reducedField);
    result(cellNumber/100, 2) = sum(abs(Te-TeSimulation)./Te) / length(reducedField);
end

plot(100:100:10000, result(:,1), 100:100:10000, result(:,2));