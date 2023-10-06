/** \file
 *
 *  Unit tests for nonuniform CAR operator
 *
 *  \author Jop Hendrikx
 *  \date   October 2023
 */

#include "LoKI-B/Grid.h"
#include "source/Operators.cpp"
#include "LoKI-B/Constant.h"
#include "LoKI-B/Gnuplot.h"
#include "LoKI-B/EedfUtilities.h"

#include "tests/TestUtilities.h"

auto json = R"json(
{
    workingConditions:
        gasPressure: 133.32                         % in Pa
        gasTemperature: 300                         % in K
        electronDensity: 1e19                       % in m-3
        electronTemperature: 1                      % in eV
        chamberLength: 1.0                          % in m
        chamberRadius: 1.0                          % in m
        reducedField: logspace(-4,2,100)            % in Td
        excitationFrequency: 0                      % in Hz

    electronKinetics:
    isOn: true             % true or false (to activate of deactivate the electron Kinetics)
    eedfType: boltzmann    % boltzmann or maxwellian
    ionizationOperatorType: usingSDCS % conservative, oneTakesAll, equalSharing or usingSDCS
    growthModelType: spatial % temporal or spatial
    includeEECollisions: false
    LXCatFiles:            % cross section files
        - Oxygen/O2_LXCat.txt               
        - Oxygen/O2_rot_LXCat.txt
    %   CARgases:             % gases for which CAR is activated
    %    - O2
    gasProperties:        % properties of the gases (S.I. Units)
        mass: Databases/masses.txt
        fraction:
        - O2 = 1
        harmonicFrequency: Databases/harmonicFrequencies.txt
        anharmonicFrequency: Databases/anharmonicFrequencies.txt
        rotationalConstant: Databases/rotationalConstants.txt
        electricQuadrupoleMoment: Databases/quadrupoleMoment.txt
        OPBParameter: Databases/OPBParameter.txt
    stateProperties:      % properties of the states (S.I. Units except for the energy [eV])
        energy:
        - O2(X,v=*) = harmonicOscillatorEnergy
        - O2(X,v=0,J=*) = rigidRotorEnergy  
        statisticalWeight:
        - O2(X,v=*) = 3.0
        - O2(X,v=0,J=*) = rotationalDegeneracy
        population:
        - O2(X) = 1.0
        - O2(X,v=0) = 1.0
    %       - O2(X,v=*) = boltzmannPopulation@gasTemperature
        - O2(X,v=0,J=*) = boltzmannPopulation@gasTemperature
    numerics: % configuration of numerical details of the simulation
        energyGrid:             % properties of the energy grid (in eV)
        maxEnergy: 1
        cellNumber: 1000
        smartGrid:            % configuration of the smart grid
            minEedfDecay: 20    % minimun number of decades of decay for the EEDF
            maxEedfDecay: 25    % maximum number of decades of decay for the EEDF
            updateFactor: 0.05  % factor used to increase or decrease the maximum value of the energy grid
        maxPowerBalanceRelError: 1e-9       % threshold for the relative power balance warning message
        nonLinearRoutines:
        algorithm: mixingDirectSolutions  % mixingDirectSolutions or iterativeSolution
        mixingParameter: 0.7              % mixingDirectSolutions mixing parameter from 0 to 1
        maxEedfRelError: 1e-9             % maximum difference for each eedf component between two iterations (stop criteria)
    %       odeSetParameters:                 % optional parameters for the ode solver of the "iterativeSolution" algorithm
    %         MaxStep: 1e-7

    % --- configuration for the heavy species kinetics ---
    chemistry:
        isOn: false

        % --- configuration of the graphical user interface ---
        gui: 
        isOn: true
        refreshFrequency: 1

    % ---  configuration of the output files ---
    output: 
    isOn: false
    folder: O2Swarm_discrete_rotations_300K
    dataFiles:
        - eedf
        - swarmParameters
        - rateCoefficients
        - powerBalance
        - lookUpTable
)json"_json;

int main()
{
    using namespace loki;

    const unsigned nCells = 100;
    const double uMax = 2; // eV
    const double T = 300; 
    const double eon = 10;
    const double won = 1;

    Vector fieldCrossSection = Vector::Ones(nCells+1);
    Vector elasticCrossSection = Vector::Ones(nCells+1);

    Grid::Vector ls(nCells+1); 
    ls << Vector::LinSpaced(nCells + 1, 0.0, 1.0);
    Grid grid1(ls,uMax,false); 
    Grid grid2(nCells, uMax);

    SparseMatrix M1(nCells,nCells);
    SparseMatrix M2(nCells,nCells);

    Vector O = Vector::Ones(10);
    CAROperator caroperator(json);
    caroperator.evaluate(grid1, T, M1);
    caroperator.evaluate(grid2, T, M2);
    test_expr(M1.isApprox(M2));

    test_report;
    return nerrors;
}
