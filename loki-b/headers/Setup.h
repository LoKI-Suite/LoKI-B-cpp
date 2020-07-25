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
     * of, and store the different sections, present in the main setup file.
     *
     * Please note that these specific classes are all defined as structs, such
     * that their member variables are directly available. This is done since
     * the only purpose of these classes will be accessing their member variables
     * to gain the required information to further parse the data.
     */

    /*
     * The SetupBase struct provides a base class that all setup structures will
     * inherit from. It contains function template 'parse' that each class
     * can override to specify how the information in this class can be obtained
     * from (a section of) the input file. The template argument is either a
     * string (legacy file format) or a JSON object. Furthermore, it defines a
     * 'parseSubStructure' template that accepts a SetupBase struct by reference,
     * which is then filled.
     */

    struct SetupBase {

	template <class SubStructure>
        static bool parseSubStructure(const std::string &content,
                const std::string &fieldName, SubStructure &subStruct);
    };

    /* ------- WORKING CONDITIONS ------- */

    /*
     * WorkingConditionsSetup is an auxiliary structure meant to store the parameters
     * from the setup file that are concerning the working conditions of the simulation.
     *
     * It is a substructure of the electronKineticsSetup class.
     */

    struct WorkingConditionsSetup : public SetupBase {
        std::string reducedField;
        std::string electronTemperature;
        double excitationFrequency = 0.;
        double gasPressure = 0.;
        double gasTemperature = 0.;
        double electronDensity = 0.;
        double chamberLength = 0.;
        double chamberRadius = 0.;

        bool parse(const std::string &sectionContent);
    };

    /* ------- GAS PROPERTIES ------- */

    /*
     * GasPropertiesSetup is an auxiliary structure that stores the properties from the setup file
     * concerning the specified gases.
     *
     * Note that its fields are all strings since gas properties are supplied through a
     * database file.
     *
     * It is a substructure of the ElectronKineticsSetup class.
     */

    struct GasPropertiesSetup : public SetupBase {
        std::string mass;
        std::vector<std::string> fraction;
        std::string harmonicFrequency,
                    anharmonicFrequency,
                    rotationalConstant,
                    electricQuadrupoleMoment,
                    OPBParameter;

        bool parse(const std::string &sectionContent);
    };

    /* ------- STATE PROPERTIES ------- */

    /*
     * StatePropertiesSetup is an auxiliary structure that stores the properties from the
     * setup file concerning the specified gases.
     *
     * Note that its fields are all strings since state properties are supplied through
     * either a file, function, or direct value. The type of the input data will be deduced
     * in the second parsing step and the parser will handle accordingly.
     *
     * It is a substructure of the ElectronKineticsSetup class.
     */

    struct StatePropertiesSetup : public SetupBase {
        std::vector<std::string> energy,
                                 statisticalWeight,
                                 population;

        bool parse(const std::string &sectionContent);
    };

    /* ------- NUMERICS ------- */

    /*
     * SmartGridSetup is a setup structure that stores the user defined settings
     * concerning the smart grid from the main setup file.
     *
     * It is a substructure of the EnergyGrid class.
     */

    struct SmartGridSetup : public SetupBase {
        uint32_t minEedfDecay = 0,
                 maxEedfDecay = 0;
        double updateFactor = 0;

        bool parse(const std::string &sectionContent);
    };

    /*
     * EnergyGridSetup is a setup structure that stores the user defined settings
     * concerning the energy grid from the main setup file.
     *
     * It is a substructure of the NumericsSetup class.
     */

    struct EnergyGridSetup : public SetupBase {
        double maxEnergy;
        uint32_t cellNumber;
        SmartGridSetup smartGrid;

        bool parse(const std::string &sectionContent);
    };

    /*
     * OdeSetParametersSetup is a setup structure that stores the user defined settings
     * concerning extra parameters needed for the iterativeSolution routine.
     *
     * It is a substructure of the NonLinearRoutinesSetup class.
     */

    struct OdeSetParametersSetup : public SetupBase {
        double maxStep;

        bool parse(const std::string &sectionContent);
    };

    /*
     * NonLinearRoutinesSetup is a setup structure that stores the user defined settings
     * concerning the non-linear algorithm and its parameters.
     *
     * It is a substructure of the NumericsSetup class.
     */

    struct NonLinearRoutinesSetup : public SetupBase {
        std::string algorithm;
        double mixingParameter,
               maxEedfRelError;
        OdeSetParametersSetup odeSetParameters;

        bool parse(const std::string &sectionContent);
    };

    /*
     * NumericsSetup is a setup structure that stores the user defined settings
     * concerning the numerical details applied in the simulation.
     *
     * It is a substructure of the ElectronKineticsSetup class.
     */

    struct NumericsSetup : public SetupBase {
        EnergyGridSetup energyGrid;
        double maxPowerBalanceRelError = -1.;
        NonLinearRoutinesSetup nonLinearRoutines;

        bool parse(const std::string &sectionContent);
    };

    /* ------- ELECTRON KINETICS ------- */

    /*
     * ElectronKineticsSetup is an auxiliary structure meant to store the parameters
     * from the setup file that are concerning the electron kinetics of the simulation.
     *
     * It is a substructure of the Setup class.
     */

    struct ElectronKineticsSetup : public SetupBase {
        bool isOn{false};
        Enumeration::EedfType eedfType;
        uint8_t shapeParameter{0};
        Enumeration::IonizationOperatorType ionizationOperatorType;
        Enumeration::GrowthModelType growthModelType;
        bool includeEECollisions{false};
        std::vector<std::string> LXCatFiles;
        std::vector<std::string> LXCatFilesExtra;
        std::vector<std::string> effectiveCrossSectionPopulations;
        std::vector<std::string> CARgases;
        GasPropertiesSetup gasProperties;
        StatePropertiesSetup stateProperties;
        NumericsSetup numerics;

        bool parse(const std::string &sectionContent);
    };

    /* ------- OUTPUT ------- */

    /*
     * OutputSetup is a setup structure meant to store the parameters that are concerning the
     * output of the simulation.
     *
     * It is a substructure of the Setup class.
     */

    struct OutputSetup : public SetupBase {
        bool isOn;
        std::string folder;
        std::vector<std::string> dataFiles;

        bool parse(const std::string &sectionContent);
    };

    /* ------- SETUP ------- */

    /*
     * The Setup structure combines the above intermediate structures to store all
     * information from the setup file in a formatted way.
     */

    class Setup : public SetupBase {

    public:

        Setup(const std::string& fname);

        WorkingConditionsSetup workingConditions;
        ElectronKineticsSetup electronKinetics;
        OutputSetup output;

        /** \todo remove this? Needed by Output when *output* is written to
         *  produce the input file. It should be possible to do this after
         *  reading without storing the result in a member.
         */
        std::string fileContent;
    private:

        bool parse(const std::string &sectionContent);
        bool parseFile(const std::string& fileName);
    };


}


#endif //LOKI_CPP_SETUP_H
