//
// Created by daan on 2-5-19.
//

#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>

#include "LoKI-B/Log.h"
#include "LoKI-B/Parse.h"
#include "LoKI-B/Setup.h"

// DONE: Think which is desirable:
//  1. Should fields that are required and cannot be parsed immediately throw an
//     exception?
//  2. Should fields that are required and cannot be parsed log a warning such
//     that their subsection (if required) can throw an exception. (Current implementation)

/*
 * The SET definition provides a shorthand for setting a member variable of a setup struct.
 * The user has to pass
 *
 * The '#property' yields the passed variable name as a string. Hence, SET requires that the
 * designated member variable has the same name as the corresponding field in the input file.
 */
#define SET(section, property)                                                                                         \
    if (!setField(section, #property, property))                                                                \
        Log<ParseFieldError>::Warning(#property);

/*
 * The R_SET variant will cause a boolean function to return false when the setField function
 * is unsuccessful. Use this variant for fields that are required to perform the simulation.
 */
#define R_SET(section, property)                                                                                       \
    if (!setField(section, #property, property))                                                                \
    {                                                                                                                  \
        Log<ParseFieldError>::Warning(#property);                                                                      \
        return false;                                                                                                  \
    }

/*
 * The SUB_STRUCT definition provides a shorthand for setting a sub structure of a SetupBase
 * object.
 */
#define SUB_STRUCT(section, subStruct) parseSubStructure(section, #subStruct, subStruct);

/*
 * The R_SUB_STRUCT requires the sub structures to be present in the input file, and parsed
 * successfully. Use this variant for sub structures that are obligatory.
 */
#define R_SUB_STRUCT(section, subStruct)                                                                               \
    if (!parseSubStructure(section, #subStruct, subStruct))                                                            \
        return false;

namespace loki
{

bool SetupBase::getFieldValue(const std::string &sectionContent, const std::string &fieldName, std::string &valueBuffer)
{
    const std::regex r(fieldName + R"(:\s*(.*[^\s\n])\s*\n*)");
    std::smatch m;

    if (!std::regex_search(sectionContent, m, r))
        return false;

    valueBuffer = m[1];

    return true;
}

bool SetupBase::getList(const std::string &sectionContent, const std::string &fieldName,
                    std::vector<std::string> &container)
{

    static const std::regex r(R"(-\s*(\S+(?: \S+)*)\s*\n*)");

    for (auto it = std::sregex_iterator(sectionContent.begin(), sectionContent.end(), r);
         it != std::sregex_iterator(); ++it)
    {

        container.emplace_back(it->str(1));
    }

    return !container.empty();
}

bool SetupBase::getSection(const std::string &fileContent, const std::string &sectionTitle, std::string &sectionBuffer)
{

    // This regular expression finds the level of a specific section. In other
    // words, it finds the number of spaces that precede the section title
    const std::regex reLevel(R"((?:^|\n)( *))" + sectionTitle);
    std::smatch m;

    if (!std::regex_search(fileContent, m, reLevel))
        return false;

    const std::string levelString = m[1];

// This regular expression matches a specific section in the input file. More accurately,
// it returns the text in between the specified section and the next section on the same
// level.
#ifdef _MSVC
    // MSVC: $ matches before \n and at the end of input
    const std::regex reSection(sectionTitle + R"(:\s*\n*([^]*?)(?:(?:\n+)" + levelString + R"(\w)|(?:$\n$)))");
#else
    // OTHER: $ matches only end of input
    const std::regex reSection(sectionTitle + R"(:\s*\n*([^]*?)(?:(?:\n+)" + levelString + R"(\w)|(?:\n*$)))");
#endif

    if (!std::regex_search(fileContent, m, reSection))
        return false;

    // Restore the original level of the content to allow for further processing.
    sectionBuffer = levelString + "  ";
    sectionBuffer += m[1];

    return true;
}


template <typename T>
bool SetupBase::setField(const std::string &sectionContent, const std::string &fieldName, T &value)
{

    std::string valueBuffer;

    if (!getFieldValue(sectionContent, fieldName, valueBuffer))
        return false;

    std::stringstream s(valueBuffer);
    return (s >> value) && s.eof();
}


/*
 * The setField function needs to behave differently when it is
 * supplied with some specific types. The types for which the
 * function needs to be specialized are:
 *
 * - bool; since these values are supplied through 'true' and
 *   'false' rather than '1' and '0'.
 * - std::string; strictly speaking, this is not necessary,
 *   but in this case we can skip the type cast since the
 *   desired type is already std::string.
 * - std::vector<std::string>; in this case the behaviour
 *   is completely different and the getList function is
 *   used.
 * - enums; whenever a new enum is defined that can occur in
 *   the input file, a new specialization needs to be written.
 *
 * Note that we could have defined a different function for every
 * type (e.g. setBool, setString, ...), however, now we can set
 * any desired value through a single function.
 */
template <>
bool SetupBase::setField<std::string>(const std::string &sectionContent, const std::string &fieldName,
                                         std::string &value)
{

    return getFieldValue(sectionContent, fieldName, value);
}

template <>
bool SetupBase::setField<bool>(const std::string &sectionContent, const std::string &fieldName, bool &value)
{
    std::string valueBuffer;

    if (!getFieldValue(sectionContent, fieldName, valueBuffer))
        return false;

    std::stringstream s(valueBuffer);
    s >> std::boolalpha >> value;

    return true;
}

template <>
bool SetupBase::setField<EedfType>(const std::string &sectionContent, const std::string &fieldName, EedfType &value)
{

    std::string valueBuffer;

    if (!getFieldValue(sectionContent, fieldName, valueBuffer))
        return false;
    value = getEedfType(valueBuffer);
    return true;
}

template <>
bool SetupBase::setField<IonizationOperatorType>(const std::string &sectionContent, const std::string &fieldName,
                                                    IonizationOperatorType &value)
{
    std::string valueBuffer;

    if (!getFieldValue(sectionContent, fieldName, valueBuffer))
        return false;
    value = getIonizationOperatorType(valueBuffer);
    return true;
}

template <>
bool SetupBase::setField<GrowthModelType>(const std::string &sectionContent, const std::string &fieldName,
                                             GrowthModelType &value)
{

    std::string valueBuffer;

    if (!getFieldValue(sectionContent, fieldName, valueBuffer))
        return false;
    value = getGrowthModelType(valueBuffer);
    return true;
}

template <>
bool SetupBase::setField<std::vector<std::string>>(const std::string &sectionContent, const std::string &fieldName,
                                                      std::vector<std::string> &value)
{
    std::string fieldContent;

    if (!getSection(sectionContent, fieldName, fieldContent))
        return false;

    return getList(fieldContent, fieldName, value);
}


/*
 * The parseSubStructure function fills a substructure with the data available in the
 * supplied (section of the) input file. This substructure is derived from the
 * SetupBase struct and it is passed by reference. Note that every such derived class
 * has overridden the 'parse' function to fill its fields and substructures as
 * desired.
 */

template <class SubStructure>
bool SetupBase::parseSubStructure(const std::string &content, const std::string &fieldName, SubStructure &subStruct)
{
    std::string sectionContent;

    /// \todo Investigate me: Added an extra new line character for MSVC regular expressions to work.
    if (getSection(content + '\n', fieldName, sectionContent))
    {
        // Same here.
        if (!subStruct.parse(sectionContent + '\n'))
        {
            Log<ParseSectionError>::Error(fieldName);
            return false;
        }
    }
    else
    {
        Log<MissingSectionError>::Warning(fieldName);
        return false;
    }

    return true;
}

Setup::Setup(const std::string &fname)
{
    if (!parseFile(fname))
    {
        throw std::runtime_error("Error parsing input file '" + fname + "'.");
    }
}

/*
 * The parseFile function is only available to the main Setup class. The user
 * passes the name of the input file. The function will then extract the text
 * from the input file, remove any comments and pass it to the parse function.
 */

bool Setup::parseFile(const std::string &fileName)
{
#ifdef EMSCRIPTEN
    const std::string inputPath{""};
#else
    const std::string inputPath("../input/");
#endif

    if (!Parse::stringBufferFromFile(inputPath + fileName, fileContent))
    {
        Log<FileError>::Error(inputPath + fileName);
        return false;
    }
    return this->parse(fileContent);
}

/*
 * Every derived class of the BaseSetup struct needs to override the parse function.
 * Inside this function the user should call setField on all the class member
 * variables, and SetupBase::parseSubStructure on all its substructures.
 *
 * At the top of this file, there are some defines that ease this process. However,
 * they do require the member variable to have the same name as specified in the
 * input file.
 *
 * Note that due to the lack of reflection/inspection in C++, we cannot reference
 * members by a string representation or loop through them. Therefore, we have to
 * call either (R_)SET() or (R_)SUB_STRUCT() for each member individually.
 */

bool Setup::parse(const std::string &sectionContent)
{

    R_SUB_STRUCT(sectionContent, workingConditions)
    R_SUB_STRUCT(sectionContent, electronKinetics)
    R_SUB_STRUCT(sectionContent, output)

    return true;
}

bool WorkingConditionsSetup::parse(const std::string &sectionContent)
{
    // TODO: Check whether 'reducedField' is present in the case that eedfType is boltzmann
    //  (and subsequently that 'electronTemperature' is present when it is prescribed).

    SET(sectionContent, reducedField)
    SET(sectionContent, electronTemperature)
    R_SET(sectionContent, excitationFrequency)
    R_SET(sectionContent, gasPressure)
    R_SET(sectionContent, gasTemperature)
    R_SET(sectionContent, electronDensity)
    R_SET(sectionContent, chamberLength)
    R_SET(sectionContent, chamberRadius)

    return true;
}

bool ElectronKineticsSetup::parse(const std::string &sectionContent)
{
    R_SET(sectionContent, isOn)
    R_SET(sectionContent, eedfType)
    SET(sectionContent, shapeParameter)
    R_SET(sectionContent, ionizationOperatorType)
    R_SET(sectionContent, growthModelType)
    R_SET(sectionContent, includeEECollisions)
    R_SET(sectionContent, LXCatFiles)
    SET(sectionContent, LXCatFilesExtra)
    SET(sectionContent, effectiveCrossSectionPopulations)
    SET(sectionContent, CARgases)

    R_SUB_STRUCT(sectionContent, gasProperties)
    R_SUB_STRUCT(sectionContent, stateProperties)
    R_SUB_STRUCT(sectionContent, numerics)

    return true;
}

bool GasPropertiesSetup::parse(const std::string &sectionContent)
{
    R_SET(sectionContent, mass)
    R_SET(sectionContent, fraction)
    R_SET(sectionContent, harmonicFrequency)
    R_SET(sectionContent, anharmonicFrequency)
    R_SET(sectionContent, rotationalConstant)
    R_SET(sectionContent, electricQuadrupoleMoment)
    R_SET(sectionContent, OPBParameter)

    return true;
}

bool StatePropertiesSetup::parse(const std::string &sectionContent)
{
    SET(sectionContent, energy)
    SET(sectionContent, statisticalWeight)
    R_SET(sectionContent, population)

    return true;
}

bool NumericsSetup::parse(const std::string &sectionContent)
{
    SET(sectionContent, maxPowerBalanceRelError)

    R_SUB_STRUCT(sectionContent, energyGrid)
    R_SUB_STRUCT(sectionContent, nonLinearRoutines)

    return true;
}

bool EnergyGridSetup::parse(const std::string &sectionContent)
{
    R_SET(sectionContent, maxEnergy)
    R_SET(sectionContent, cellNumber)

    SUB_STRUCT(sectionContent, smartGrid)

    return true;
}

bool SmartGridSetup::parse(const std::string &sectionContent)
{
    R_SET(sectionContent, minEedfDecay)
    R_SET(sectionContent, maxEedfDecay)
    R_SET(sectionContent, updateFactor)

    return true;
}

bool NonLinearRoutinesSetup::parse(const std::string &sectionContent)
{
    R_SET(sectionContent, algorithm)
    R_SET(sectionContent, mixingParameter)
    R_SET(sectionContent, maxEedfRelError)

    SUB_STRUCT(sectionContent, odeSetParameters)

    return true;
}

bool OdeSetParametersSetup::parse(const std::string &sectionContent)
{
    R_SET(sectionContent, maxStep)

    return true;
}

bool OutputSetup::parse(const std::string &sectionContent)
{
    R_SET(sectionContent, isOn)
    R_SET(sectionContent, folder)
    R_SET(sectionContent, dataFiles)

    return true;
}
} // namespace loki
