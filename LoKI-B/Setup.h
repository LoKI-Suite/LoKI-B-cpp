//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_SETUP_H
#define LOKI_CPP_SETUP_H

#include <string>
#include <vector>

#include "LoKI-B/Enumeration.h"
#include "LoKI-B/Exports.h"

namespace loki
{
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

struct lokib_export SetupBase
{

    template <class SubStructure>
    static bool parseSubStructure(const std::string &content, const std::string &fieldName, SubStructure &subStruct);
protected:
    /** The getFieldValue function will extract the value of a given field name.
     *  From a section in the input file. This function is specifically designed
     *  to extract single line values (thus they are not a section or list). The
     *  string buffer to hold the value is passed by reference and the function
     *  returns a boolean to specify whether the operation was successful.
     *
     *  NOTE: This function does not deal with ambiguous field names. E.g. the "isOn"
     *  field occurs multiple times in the file. The user is advised to only retrieve
     *  fields that are one level above the level of "sectionContent". E.g. retrieve
     *  "isOn" when "sectionContent" contains the contents of the "electronKinetics"
     *  section.
     */
    bool getFieldValue(const std::string &sectionContent, const std::string &fieldName, std::string &valueBuffer);
    /** The getList function retrieves all entries in a list type field, e.g.:
     *
     *  dataFiles:
     *    - eedf
     *    - swarmParameters
     *
     *  and returns them as a vector of strings (thus {"eedf", "swarmParameters"}
     *  in the example). This vector is passed by reference as an argument. and
     *  the function returns a boolean to specify whether the operation was
     *  successful.
     */
    bool getList(const std::string &sectionContent, const std::string &fieldName,
                        std::vector<std::string> &container);
    /** getSection is a function that retrieves the contents of a specified section
     *  and stores them in the "sectionBuffer" string. Furthermore, it returns a boolean
     *  to indicate whether the operation was successful.
     */
    static bool getSection(const std::string &fileContent, const std::string &sectionTitle, std::string &sectionBuffer);

    /** The setField function extracts a value from a field in the input file,
     *  casts it to the appropriate type and assigns it to the 'value' argument
     *  which is passed by reference. It returns a boolean to give an indication
     *  whether the operation was successful.
     *
     *  Note that this is a template function, and the type of 'value' is arbitrary,
     *  a standard string can be converted to most basic types using a stringstream.
     *  However, for some basic types (e.g. bool) and custom types, such as enums,
     *  we can specialize this function to alter its behaviour.
     *
     *  The definitions of these specializations can be found below the Parse struct.
     */
    template <typename T>
    bool setField(const std::string &sectionContent, const std::string &fieldName, T &value);
};

/* ------- WORKING CONDITIONS ------- */

/*
 * WorkingConditionsSetup is an auxiliary structure meant to store the parameters
 * from the setup file that are concerning the working conditions of the simulation.
 *
 * It is a substructure of the electronKineticsSetup class.
 */

struct lokib_export WorkingConditionsSetup : public SetupBase
{
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

struct GasPropertiesSetup : public SetupBase
{
    std::string mass;
    std::vector<std::string> fraction;
    std::string harmonicFrequency, anharmonicFrequency, rotationalConstant, electricQuadrupoleMoment, OPBParameter;

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

struct lokib_export StatePropertiesSetup : public SetupBase
{
    std::vector<std::string> energy, statisticalWeight, population;

    bool parse(const std::string &sectionContent);
};

/* ------- NUMERICS ------- */

/*
 * SmartGridSetup is a setup structure that stores the user defined settings
 * concerning the smart grid from the main setup file.
 *
 * It is a substructure of the EnergyGrid class.
 */

struct lokib_export SmartGridSetup : public SetupBase
{
    uint32_t minEedfDecay = 0, maxEedfDecay = 0;
    double updateFactor = 0;

    bool parse(const std::string &sectionContent);
};

/*
 * EnergyGridSetup is a setup structure that stores the user defined settings
 * concerning the energy grid from the main setup file.
 *
 * It is a substructure of the NumericsSetup class.
 */

struct lokib_export EnergyGridSetup : public SetupBase
{
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

struct lokib_export OdeSetParametersSetup : public SetupBase
{
    double maxStep;

    bool parse(const std::string &sectionContent);
};

/*
 * NonLinearRoutinesSetup is a setup structure that stores the user defined settings
 * concerning the non-linear algorithm and its parameters.
 *
 * It is a substructure of the NumericsSetup class.
 */

struct lokib_export NonLinearRoutinesSetup : public SetupBase
{
    std::string algorithm;
    double mixingParameter, maxEedfRelError;
    OdeSetParametersSetup odeSetParameters;

    bool parse(const std::string &sectionContent);
};

/*
 * NumericsSetup is a setup structure that stores the user defined settings
 * concerning the numerical details applied in the simulation.
 *
 * It is a substructure of the ElectronKineticsSetup class.
 */

struct lokib_export NumericsSetup : public SetupBase
{
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

struct lokib_export ElectronKineticsSetup : public SetupBase
{
    bool isOn{false};
    EedfType eedfType;
    uint8_t shapeParameter{0};
    IonizationOperatorType ionizationOperatorType;
    GrowthModelType growthModelType;
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

struct lokib_export OutputSetup : public SetupBase
{
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

class lokib_export Setup : public SetupBase
{

  public:
    Setup(const std::string &fname);

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
    bool parseFile(const std::string &fileName);
};

} // namespace loki

#endif // LOKI_CPP_SETUP_H
