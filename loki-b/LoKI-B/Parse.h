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
    std::stringstream ss{valueString};
    return static_cast<bool>(ss >> value) && ss.eof();
}

/** Parse \a valueString into a double and return the result. If the conversion
 *  fails or is incomplete (trailing characters are present), a
 *  std::runtime_error is thrown.
 */
inline double getValue(const std::string &valueString)
{
    double value;
    if (getValue(valueString,value))
    {
        return value;
    }
    else
    {
        throw std::runtime_error("Error converting string '" + valueString + "' to a number.");
    }
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

} // namespace Parse
} // namespace loki

#endif // LOKI_CPP_PARSE_H
