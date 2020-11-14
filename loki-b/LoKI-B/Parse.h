//
// Created by daan on 6-5-19.
//

#ifndef LOKI_CPP_PARSE_H
#define LOKI_CPP_PARSE_H

#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>
#include <string>

#include "LoKI-B/Enumeration.h"
#include "LoKI-B/StandardPaths.h"
#include "LoKI-B/json.h"

namespace loki {
namespace Parse {

    inline bool isNumerical(const std::string &str)
    {
        static const std::regex reNum(R"(\s*\d*\.?\d+\s*)");

        return std::regex_match(str, reNum);
    }

    /** Parse \a valueString into double \a value. The boolean return
     *  value returns tru if the conversion was successful, false otherwise.
     */
    inline bool getValue(const std::string &valueString, double &value)
    {
        //            const std::regex r(R"(\s*(\d+\.?\d*)\s*\n*)");
        //            std::smatch m;
        //
        //            if (!std::regex_match(valueString, r)) return false;

        std::stringstream ss(valueString);

        return static_cast<bool>(ss >> value);
    }

    /*
     * The getFieldValue function will extract the value of a given field name.
     * From a section in the input file. This function is specifically designed
     * to extract single line values (thus they are not a section or list). The
     * string buffer to hold the value is passed by reference and the function
     * returns a boolean to specify whether the operation was successful.
     *
     * NOTE: This function does not deal with ambiguous field names. E.g. the "isOn"
     * field occurs multiple times in the file. The user is advised to only retrieve
     * fields that are one level above the level of "sectionContent". E.g. retrieve
     * "isOn" when "sectionContent" contains the contents of the "electronKinetics"
     * section.
     */
    inline bool getFieldValue(const std::string &sectionContent, const std::string &fieldName, std::string &valueBuffer)
    {
        const std::regex r(fieldName + R"(:\s*(.*[^\s\n])\s*\n*)");
        std::smatch m;

        if (!std::regex_search(sectionContent, m, r))
            return false;

        valueBuffer = m[1];

        return true;
    }

    /*
     * The getList function retrieves all entries in a list type field, e.g.:
     *
     * dataFiles:
     *   - eedf
     *   - swarmParameters
     *
     * and returns them as a vector of strings (thus {"eedf", "swarmParameters"}
     * in the example). This vector is passed by reference as an argument. and
     * the function returns a boolean to specify whether the operation was
     * successful.
     */
    inline bool getList(const std::string &sectionContent, const std::string &fieldName,
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

    /*
     * getSection is a function that retrieves the contents of a specified section
     * and stores them in the "sectionBuffer" string. Furthermore, it returns a boolean
     * to indicate whether the operation was successful.
     */
    inline bool getSection(const std::string &fileContent, const std::string &sectionTitle, std::string &sectionBuffer)
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

    /*
     * removeComments is a function that takes a string as an argument. It strips the
     * string from comments and returns it.
     */
    inline std::string removeComments(const std::string &content)
    {
        std::string content_clean;

        static const std::regex reLine(R"(\n\s*%[^\n]*)");
        content_clean = std::regex_replace(content, reLine, "");

        static const std::regex reClean(R"(%[^]*?(?:\n|$))");
        return std::regex_replace(content_clean, reClean, "\n");
    }

    /* -- statePropertyDataType --
     * Deduces whether an entry in the stateProperties section describes loading
     * of state properties by direct value, file or function. The result is
     * returned as a StatePropertyDataType enumeration type.
     */

    inline StatePropertyDataType statePropertyDataType(const std::string &propertyString, std::string &buffer)
    {
        static const std::regex reProperty(R"(.*=\s*(.+?)\s*$)");
        static const std::regex reValue(R"([\.0-9]+)");

        std::smatch m;

        if (std::regex_search(propertyString, m, reProperty))
        { // value or function
            buffer = m.str(1);

            if (std::regex_match(buffer, reValue))
            {
                return StatePropertyDataType::direct;
            }

            return StatePropertyDataType::function;
        }

        buffer = propertyString;

        return StatePropertyDataType::file;
    }

    /* -- propertyFunctionAndArguments --
     * Extracts the property function name and arguments from a given string and
     * stores them in two separate strings.
     */

    inline bool propertyFunctionAndArguments(const std::string &totalString, std::string &functionName,
                                             std::string &argumentString)
    {
        static const std::regex reFuncArgs(R"(\s*(\w+)@?(.*))");
        std::smatch m;

        if (!std::regex_match(totalString, m, reFuncArgs))
            return false;

        functionName = m.str(1);
        argumentString = m.str(2);

        return true;
    }

    /* -- argumentsFromString --
     * Extracts the separate arguments from a string containing the arguments, finds
     * their corresponding double values and pushes them into the arguments vector. E.g.
     * "[gasTemperature, 1000]" first yields "gasTemperature" which is looked up in
     * the argumentMap to obtain its value, and then yields "1000" which is converted
     * into a double and pushed into the arguments vector.
     */

    inline bool argumentsFromString(const std::string &argumentString, std::vector<double> &arguments,
                                    const std::map<std::string, double *> &argumentMap)
    {
        static const std::regex r(R"(\s*([\w\.]+)\s*(?:[,\]]|$))");

        for (auto it = std::sregex_iterator(argumentString.begin(), argumentString.end(), r);
             it != std::sregex_iterator(); ++it)
        {
            std::string current = it->str(1);

            if (isNumerical(current))
            {
                double value;

                if (!getValue(current, value))
                    return false;

                arguments.emplace_back(value);
            }
            else
            {
                if (argumentMap.count(current) == 0)
                    return false;

                arguments.emplace_back(*argumentMap.at(current));
            }
        }

        return true;
    }

    /* -- stringBufferFromFile --
     * Loads the complete content of a specified file into the given std::string.
     */

    inline bool stringBufferFromFile(const std::string &fileName, std::string &buffer)
    {
        std::ifstream in(INPUT "/" + fileName);

        if (!in)
            return false;

        std::stringstream ss;
        ss << in.rdbuf();

        buffer = removeComments(ss.str());

        return true;
    }

    /* -- gasProperty --
     * Tries to parse the value of a given property for a given gas. The 'content' string
     * should store the content of the database file corresponding to the passed property
     * (i.e. mass or anharmonicFrequency).
     */

    inline bool gasProperty(const std::string &gasName, double &property, const std::string &content)
    {
        const std::regex r(R"((?:^|\n))" + gasName + R"(\s+(.*)\s*)");
        std::smatch m;

        if (!std::regex_search(content, m, r))
            return false;

        std::stringstream ss(m[1]);
        return static_cast<bool>(ss >> property);
    }

    /*
     * The setField function extracts a value from a field in the input file,
     * casts it to the appropriate type and assigns it to the 'value' argument
     * which is passed by reference. It returns a boolean to give an indication
     * whether the operation was successful.
     *
     * Note that this is a template function, and the type of 'value' is arbitrary,
     * a standard string can be converted to most basic types using a stringstream.
     * However, for some basic types (e.g. bool) and custom types, such as enums,
     * we can specialize this function to alter its behaviour.
     *
     * The definitions of these specializations can be found below the Parse struct.
     */

    template <typename T>
    bool setField(const std::string &sectionContent, const std::string &fieldName, T &value)
    {

        std::string valueBuffer;

        if (!getFieldValue(sectionContent, fieldName, valueBuffer))
            return false;

        std::stringstream s(valueBuffer);
        s >> value;

        return true;
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
 *
 * Use "inline" with fully specialized templates to refrain from
 * breaking the "One Definition Rule". An alternative is to add
 * a Parse.cpp file and specify the specializations there.
 */
template <>
inline bool setField<std::string>(const std::string &sectionContent, const std::string &fieldName,
                                         std::string &value)
{

    return getFieldValue(sectionContent, fieldName, value);
}

template <>
inline bool setField<bool>(const std::string &sectionContent, const std::string &fieldName, bool &value)
{
    std::string valueBuffer;

    if (!getFieldValue(sectionContent, fieldName, valueBuffer))
        return false;

    std::stringstream s(valueBuffer);
    s >> std::boolalpha >> value;

    return true;
}

template <>
inline bool setField<EedfType>(const std::string &sectionContent, const std::string &fieldName, EedfType &value)
{

    std::string valueBuffer;

    if (!getFieldValue(sectionContent, fieldName, valueBuffer))
        return false;
    value = getEedfType(valueBuffer);
    return true;
}

template <>
inline bool setField<IonizationOperatorType>(const std::string &sectionContent, const std::string &fieldName,
                                                    IonizationOperatorType &value)
{
    std::string valueBuffer;

    if (!getFieldValue(sectionContent, fieldName, valueBuffer))
        return false;
    value = getIonizationOperatorType(valueBuffer);
    return true;
}

template <>
inline bool setField<GrowthModelType>(const std::string &sectionContent, const std::string &fieldName,
                                             GrowthModelType &value)
{

    std::string valueBuffer;

    if (!getFieldValue(sectionContent, fieldName, valueBuffer))
        return false;
    value = getGrowthModelType(valueBuffer);
    return true;
}

template <>
inline bool setField<std::vector<std::string>>(const std::string &sectionContent, const std::string &fieldName,
                                                      std::vector<std::string> &value)
{
    std::string fieldContent;

    if (!getSection(sectionContent, fieldName, fieldContent))
        return false;

    return getList(fieldContent, fieldName, value);
}

} // namespace Parse
} // namespace loki

#endif // LOKI_CPP_PARSE_H
