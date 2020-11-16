//
// Created by daan on 6-5-19.
//

#ifndef LOKI_CPP_PARSE_H
#define LOKI_CPP_PARSE_H

#include <string>
#include <sstream>
#include <fstream>

#include "LoKI-B/StandardPaths.h"

namespace loki {
namespace Parse {

/** Parse \a valueString into double \a value. The boolean return
 *  value returns tru if the conversion was successful, false otherwise.
 */
inline bool getValue(const std::string &valueString, double &value)
{
    std::stringstream ss{valueString};
    return (ss >> value) && ss.eof();
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

inline std::string& searchAndReplaceInPlace(
                                std::string& str,
                                const std::string& key,
                                const std::string& repl)
{
        const std::string::size_type sz = key.length();
        if (!sz)
        {
                throw std::runtime_error("search and replace: empty search string provided.");
        }
        for (std::string::size_type pos=0;                    // start looking at position 0
             (pos = str.find(key,pos)) != std::string::npos;  // update pos, stop if there is no more match
             pos += repl.length() )                           // advance to the end of the replaced substring
        {
                str.replace(pos, sz, repl);
        }
        return str;
}

inline std::string searchAndReplaceCopy(
                                const std::string& str,
                                const std::string& key,
                                const std::string& repl)
{
        std::string copy(str);
        return searchAndReplaceInPlace(copy,key,repl);
}

/** Reads characters from stream \a is into \a buffer.
 *  Comments and empty lines are removed.
 */
inline bool removeComments(std::istream& is, std::string &buffer)
{
    std::string line;
    while (!is.eof())
    {
        std::getline(is,line);
        // remove everything from % onward (if one is found).
        std::string::size_type comment_pos = line.find('%');
        line = line.substr(0,comment_pos);
        // remove trailing whitespace
        line.erase(line.find_last_not_of(" \n\r\t")+1);
        // If there is anything left, append it to the buffer.
        // Prepend a newline character if buffer already contains
        // something.
        if (!line.empty())
        {
            if (!buffer.empty())
            {
                buffer += '\n';
            }
            buffer += line;
        }
    }
    return true;
}

inline bool stringBufferFromFile(const std::string &fileName, std::string &buffer)
{
    std::ifstream is(INPUT "/" + fileName);
    if (!is)
    {
        return false;
    }
    return removeComments(is,buffer);
}

} // namespace Parse
} // namespace loki

#endif // LOKI_CPP_PARSE_H
