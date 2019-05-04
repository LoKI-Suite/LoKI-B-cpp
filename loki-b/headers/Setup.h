//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_SETUP_H
#define LOKI_CPP_SETUP_H

#include <string>
#include <vector>
#include <Enumeration.h>

namespace loki {
    /*
     * All structures presented in this header file are used to do a first parse
     * of, and store the different sections in the main setup file.
     *
     * Please note that these specific classes are all defined as structs, such
     * that their member variables are directly available. This is done since
     * the only purpose of these classes will be accessing their member variables
     * to gain the required information to further parse the data.
     */

    /* ------- WORKING CONDITIONS ------- */

    /*
     * WorkingConditionsSetup is an auxiliary structure meant to store the parameters
     * from the setup file that are concerning the working conditions of the simulation.
     *
     * It is a substructure of the electronKineticsSetup class.
     */

    struct WorkingConditionsSetup {
        std::string reducedField;
        std::string electronTemperature;
        double excitationFrequency = 0.;
        double gasPressure = 0.;
        double gasTemperature = 0.;
        double electronDensity = 0.;
        double chamberLength = 0.;
        double chamberRadius = 0.;
    };

    /* ------- GAS PROPERTIES ------- */

    /*
     * GasPropertiesSetup is an auxiliary structure that stores the properties from the setup file
     * concerning the specified gasses.
     *
     * Note that its fields are all strings since gas properties are supplied through a
     * database file.
     *
     * It is a substructure of the ElectronKineticsSetup class.
     */

    struct GasPropertiesSetup {
        std::string mass;
        std::vector<std::string> fraction;
        std::string harmonicFrequency,
                    anharmonicFrequency,
                    rotationalConstant,
                    electricQuadrupoleMoment,
                    OPBParameter;

    };

    /* ------- STATE PROPERTIES ------- */

    /*
     * StatePropertiesSetup is an auxiliary structure that stores the properties from the
     * setup file concerning the specified gasses.
     *
     * Note that its fields are all strings since gas properties are supplied through either
     * a file, function, or direct value.
     *
     * It is a substructure of the ElectronKineticsSetup class.
     */

    // TODO: How do we want to store the state entries? Since they can have different types
    //  --> (?) Create a simple base class that stores the type and has methods to access the data
    //  inherit this class into three subclasses corresponding to a file, function and direct value.

    struct StatePropertiesSetup {


    };

    /* ------- NUMERICS ------- */

    /*
     * SmartGridSetup is a setup structure that stores the user defined settings
     * concerning the smart grid from the main setup file.
     *
     * It is a substructure of the EnergyGrid class.
     */

    struct SmartGridSetup {
        uint32_t minEedfDecay,
                 maxEedfDecay;
        double updateFactor;
    };

    /*
     * EnergyGridSetup is a setup structure that stores the user defined settings
     * concerning the energy grid from the main setup file.
     *
     * It is a substructure of the NumericsSetup class.
     */

    struct EnergyGridSetup {
        double maxEnergy;
        uint32_t cellNumber;
        SmartGridSetup smartGrid;
    };

    /*
     * OdeSetParametersSetup is a setup structure that stores the user defined settings
     * concerning extra parameters needed for the iterativeSolution routine.
     *
     * It is a substructure of the NonLinearRoutinesSetup class.
     */

    struct OdeSetParametersSetup {
        double maxStep;
    };

    /*
     * NonLinearRoutinesSetup is a setup structure that stores the user defined settings
     * concerning the non-linear algorithm and its parameters.
     *
     * It is a substructure of the NumericsSetup class.
     */

    struct NonLinearRoutinesSetup {
        std::string algorithm;
        double mixingParameter,
               maxEedfRelError;
        OdeSetParametersSetup odeSetParameters;
    };

    /*
     * NumericsSetup is a setup structure that stores the user defined settings
     * concerning the numerical details applied in the simulation.
     *
     * It is a substructure of the ElectronKineticsSetup class.
     */

    struct NumericsSetup {
        EnergyGridSetup energyGrid;
        NonLinearRoutinesSetup nonLinearRoutines;
    };

    /* ------- ELECTRON KINETICS ------- */

    /*
     * ElectronKineticsSetup is an auxiliary structure meant to store the parameters
     * from the setup file that are concerning the electron kinetics of the simulation.
     *
     * It is a substructure of the Setup class.
     */

    struct ElectronKineticsSetup {
        bool isEnabled = false;
        Enumeration::EedfType eedf;
        uint8_t shapeParameter = 0;
        Enumeration::IonizationOperatorType ionizationOperator;
        Enumeration::GrowthModelType growthModel;
        bool includeEECollisions = false;
        std::vector<std::string> LXCatFiles;
        std::vector<std::string> extraLXCatFiles;
        std::vector<std::string> effectiveCrossSectionPopulations;
        std::vector<std::string> CARgases;
        GasPropertiesSetup gasProperties;
        StatePropertiesSetup stateProperties;
        NumericsSetup numerics;
    };

    /* ------- SETUP ------- */

    /*
     * The Setup structure combines the above intermediate structures to store all
     * information from the setup file in a formatted way.
     */

    struct Setup {
        WorkingConditionsSetup workingConditions;
        ElectronKineticsSetup electronKinetics;

        explicit Setup(const std::string& fileName);
    };


}


#endif //LOKI_CPP_SETUP_H
