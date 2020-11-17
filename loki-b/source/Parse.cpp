/** \file
 *
 *  Implementation of string parsing and utility functions for LoKI-B.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2020 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
 *  M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
 *  L. L. Alves
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  \author Daan Boer and Jan van Dijk (C++ version)
 */

#include "LoKI-B/Parse.h"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>

/// \todo Remove me once INPUT is no longer used.
#include "LoKI-B/StandardPaths.h"

namespace loki {
namespace Parse {

/** Parse \a valueString into double \a value. The boolean return
 *  value returns tru if the conversion was successful, false otherwise.
 */
bool getValue(const std::string &valueString, double &value)
{
    std::stringstream ss{valueString};
    return (ss >> value) && ss.eof();
}

/** Parse \a valueString into a double and return the result. If the conversion
 *  fails or is incomplete (trailing characters are present), a
 *  std::runtime_error is thrown.
 */
double getValue(const std::string &valueString)
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

/** Replace all occurrences of \a key in \a str with \a repl.
 *  The function correctly handles the case that \a key is a substring of \a repl.
 *  As an example, if in the string "AA" all "A" are to be replaced with "AB, the
 *  result will be "ABAB": the 'produced' characters 'A' will not be replaced
 *  recursively (which would result in an endless loop).
 *
 *  If the search string \a key is empty, a std::runtime_error is thrown.
 *
 *  \sa     searchAndReplaceCopy
 *  \author Jan van Dijk
 *  \date   May 2013
 */
std::string& searchAndReplaceInPlace(
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

/** Replace all occurrences of \a key in a copy of \a str with \a repl
 *  and return the result.
 *
 *  This function is similar to search_and_replace, but accepts a constant
 *  string reference. The function makes a copy of \a str, then calls
 *  searchAndReplaceInPlace to do the substitutions in that copy,
 *  and returns the result.
 *
 *  \sa     searchAndReplaceInPlace
 *  \author Jan van Dijk
 *  \date   January 2014
 */
std::string searchAndReplaceCopy(
                            const std::string& str,
                            const std::string& key,
                            const std::string& repl)
{
        std::string copy(str);
        return searchAndReplaceInPlace(copy,key,repl);
}

/** Reads characters from stream \a is into \a dest.
 *  Om entry, the result \a dest is cleared. Then the function starts
 *  extracts characters from \a is line by line, using the std::getine.
 *  For every line, it scans for a '%' character, and if one is found it
 *  removes everything starting from that point. Next, trailing whitespace
 *  is removed. When the resulting line is not empty, it is appended to the
 *  \a result (A newline character is prepanded if \a dest is not empty).
 *
 *  The function results a boolean that indicates success. At present, true
 *  is always returned, but a future extension could do additional syntax
 *  checking and return false if an error is detected. (Throwing an exception
 *  with a clear error message will be more useful in that case, though.)
 *
 *  \todo Should a % that immediately follows a non-whitespace
 *        character start a comment? Or should it be kept in a
 *        situation like "fraction: 100%"?
 *
 *  \author Daan Boer and Jan van Dijk
 *  \date   November 2020
 */
bool removeComments(std::istream& is, std::string &dest)
{
    dest.clear();
    std::string line;
    while (!is.eof())
    {
        std::getline(is,line);
        std::string::size_type comment_pos = line.find('%');
        line = line.substr(0,comment_pos);
        line.erase(line.find_last_not_of(" \n\r\t")+1);
        if (!line.empty())
        {
            if (!dest.empty())
            {
                dest += '\n';
            }
            dest += line;
        }
    }
    return true;
}

/** Open file \a filename for reading, call removeComments() to remove
 *  comments, trailing whitespace and empty lines from the stream, and
 *  send the output to \a dest.
 *
 *  The function results a boolean that indicates success. When the file
 *  cannot be opened, false is returned, otherwise the result of the call
 *  to removeComments is returned.
 *
 *  \todo the string INPUT "/" is prepended unconditionally. This is better
 *        decided by the caller. It would be cleaner and more general to
 *        not do such magic and use \a fileName as-is. The comments above
 *        already refelect that future situation.
 *
 *  \author Daan Boer and Jan van Dijk
 *  \date   November 2020
 */
bool stringBufferFromFile(const std::string &fileName, std::string &dest)
{
    std::ifstream is(INPUT "/" + fileName);
    return is && removeComments(is,dest);
}

} // namespace Parse
} // namespace loki
