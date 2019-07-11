path = '../../build/Output/delta_inelastic/eedf.txt';

fid = fopen(path, 'r');
fgetl(fid);

data = fscanf(fid,'%e',[3 inf]);

plot (data(1,:), data(2,:));