/** \file
 *
 *  Declarations of string parsing and utility functions for LoKI-B.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2024 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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
 *  \date   May 2019
 */

#ifndef LOKI_CPP_PARSE_H
#define LOKI_CPP_PARSE_H

#include "LoKI-B/Exports.h"
#include "LoKI-B/json.h"
#include <string>
#include <iosfwd>
#include <filesystem>

namespace loki {
namespace Parse {

/** An object of type FunctionCall parses a function with a particular arity
 *  from a string. The string must be of the form
 *  name '(' [ argument (',' argument )*] ')'. The string 'name' must have
 *  the same format as a function name in C++ or MATLAB, 'argument' is either
 *  a 'name' or a numerical value. Upon construction, the function name and
 *  the arguments are stored in data members for later access. The arguments
 *  are converted to JSON values; if such conversion fails it is stored as
 *  a JSON string.
 *
 *  Examples:
 *
 *    logspan(-3.0,3.0,7)
 *    _g12(x,false)
 *
 *  \sa     makeJsonFromFunctionCall
 *  \author Jan van Dijk
 *  \date   April 2024
 */
struct FunctionCall
{
    FunctionCall(unsigned arity, const std::string& str);
    using Arguments = std::vector<json_type>;
    unsigned arity() const { return m_args.size(); }
    const std::string& name() const { return m_name; }
    const Arguments& args() const { return m_args; }
    const json_type& arg(unsigned i) const { return m_args[i]; }
private:
    std::string m_name;
    Arguments m_args;
};

/** Create a JSON object from a function call string and return the result.
 *  The arguments \a fkey and \a akeys are the keys that are used for the
 *  name of the function, as parsed from \a str, and its arguments.
 *
 *  As an example, when this function is called with arguments "atan2(4,3)",
 *  "function" and {"y","x}, it will return the JSON object
 *  { "function": "atan2", "x": "3", "y": "4" }.
 *
 *  \sa     FunctionCall
 *  \author Jan van Dijk
 *  \date   April 2024
 */
json_type makeJsonFromFunctionCall(
    const std::string& str,
    const std::string& fkey,
    const std::vector<std::string>& akeys);

/** Parse \a valueString into double \a value. The function returns true if
 *  the conversion was successful, false otherwise. In the latter case the
 *  \a value is undefined after returning from this function.
 */
lokib_export bool getValue(const std::string &valueString, double &value);

/** Parse \a valueString into a double and return the result. If the conversion
 *  fails or is incomplete (trailing characters are present), a
 *  std::runtime_error is thrown.
 */
lokib_export double getValue(const std::string &valueString);

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
lokib_export std::string& searchAndReplaceInPlace(
                            std::string& str,
                            const std::string& key,
                            const std::string& repl);

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
lokib_export std::string searchAndReplaceCopy(
                            const std::string& str,
                            const std::string& key,
                            const std::string& repl);

/** Reads characters from stream \a is into \a dest.
 *  Om entry, the result \a dest is cleared. Then the function starts
 *  extracting characters from \a is line by line, using std::getline.
 *  For every line, it scans for a '%' character, and if one is found it
 *  removes everything starting from that point. Next, trailing whitespace
 *  is removed. When the resulting line is not empty, it is appended to the
 *  \a result (In this case, a newline character is prepended if \a dest is
 *  not empty --- we already  a line before).
 *
 *  The function returns a boolean that indicates success. At present, true
 *  is always returned, but a future extension could do additional syntax
 *  checking and return false if an error is detected. (Throwing an exception
 *  with a clear error message will be more useful in that case, though.)
 *
 *  Note that this function follows the MATLAB version of Loki-B in interpreting
 *  a % character as the start of a comment even if it is part of a string or if
 *  it follows a non-whitespace character.
 *
 *  \author Daan Boer and Jan van Dijk
 *  \date   November 2020
 */
lokib_export bool removeComments(std::istream& is, std::string &dest);

/** Open file \a fileName for reading, call removeComments() to remove
 *  comments, trailing whitespace and empty lines from the stream, and
 *  write the output into \a dest.
 *
 *  The function results a boolean that indicates success. When the file
 *  cannot be opened, false is returned, otherwise the result of the call
 *  to removeComments is returned.
 *
 *  \author Daan Boer and Jan van Dijk
 *  \date   November 2020
 */
lokib_export bool stringBufferFromFile(const std::filesystem::path &fileName, std::string &dest);

} // namespace Parse
} // namespace loki

#endif // LOKI_CPP_PARSE_H
