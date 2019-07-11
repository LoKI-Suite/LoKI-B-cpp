kb = 1.38064852e-23;           % Boltzmann constant in J/K
e = 1.6021766208e-19;    % Electron charge in C
me = 9.10938356e-31;        % Electron mass in Kg
unifiedAtomicMass = 1.660539040e-27;  % Unified Atomic Mass unit (UAM) in kg
bohrRadius = 5.2917721067e-11;        % Bohr radius in m
vacuumPermittivity = 8.854187817e-12; % Vacuum permittivity in F/m
plank = 6.626070040e-34;              % Plank constant in J s
speedOfLight = 299792458;             % Speed of light in vacuum in m/s
atmosphereInPa = 101325;              % Standard atmosphere in Pa
atmosphereInTorr = 760;               % Standard atmosphere in Torr (not obtained from NIST database)
energies = linspace(0, 30, 3000);

crossSection = sqrt(me ./ (2 * e .* energies)) .* 1e-12;
crossSection(1) = 0;

file = fopen('rawElastic.txt', 'w');

for i = 1:length(energies)
    fprintf(file, '%.16e\t%.16e\n', energies(i), crossSection(i));
end

fclose(file);

reducedFieldSI = 10e-21;
n = 133.32 / (kb * 300);
fieldSI = reducedFieldSI * n;

Te = 300 + e / (3 * kb * fieldSI) * (reducedFieldSI / 1e-12)^2;

energiesSI = energies .* e;

maxwellian = 2 * (kb * Te)^-3/2 .* sqrt(energiesSI ./ pi) .* exp(-energiesSI ./ (kb * Te));

semilogy(energies, maxwellian);