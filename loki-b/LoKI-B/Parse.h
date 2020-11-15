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
    std::stringstream ss(valueString);
    return static_cast<bool>(ss >> value) && ss.eof();
}

/** Remove empty lines, comments and trailing whitespace from
 *  the character stream \a content and return the result.
 */
inline std::string removeComments(const std::string &content)
{
    static const std::regex reLine(R"(\n\s*%[^\n]*)");
    static const std::regex reClean(R"(%[^]*?(?:\n|$))");
    // next line: also remove whitespace before a %
    //static const std::regex reClean(R"(\s*%[^]*?(?:\n|$))");
    std::string result;
    result = std::regex_replace(content, reLine, "");
    result = std::regex_replace(result, reClean, "\n");
    return result;
}

/** Loads the complete content of a specified file into the given std::string.
 *  Comments are removed.
 */
inline bool stringBufferFromFile(const std::string &fileName, std::string &buffer)
{
    std::ifstream is(INPUT "/" + fileName);
    if (!is)
    {
        return false;
    }
    const std::string str{
                std::istreambuf_iterator<char>(is),
                std::istreambuf_iterator<char>()};
    buffer = removeComments(str);
    return true;
}

/** Deduce whether an entry in the stateProperties section describes loading
 *  of state properties by direct value, file or function. The result is
 *  returned as a StatePropertyDataType enumeration type.
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

/** Extracts the property function name and arguments from a given string and
 *  stores them in two separate strings.
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

/** Extracts the separate arguments from a string containing the arguments, finds
 *  their corresponding double values and pushes them into the arguments vector. E.g.
 *  "[gasTemperature, 1000]" first yields "gasTemperature" which is looked up in
 *  the argumentMap to obtain its value, and then yields "1000" which is converted
 *  into a double and pushed into the arguments vector.
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

} // namespace Parse
} // namespace loki

#endif // LOKI_CPP_PARSE_H
