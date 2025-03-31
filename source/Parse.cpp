/** \file
 *
 *  Implementation of string parsing and utility functions for LoKI-B.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2025 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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

#include "LoKI-B/Parse.h"
#include <string>
#include <iostream>
#include <regex>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace loki {
namespace Parse {

FunctionCall::FunctionCall(unsigned arity, const std::string& str)
{
    /* 1. Compile a regular expression for a function with 'arity' arguments.
     */
    /* \todo Move these also outside of this function, so they
     *        can be reused in other contexts, if necessary: in
     *        particular the double-parser is non-trivial.
     */
    const std::string ws = "\\s*";
    const std::string name = "[A-Za-z_][A-Za-z_0-9]*";
    const std::string flt = "(?:[-+]?)(?:[0-9]*.)?[0-9]+(?:[eE][-+]?\\d+)?";
    const std::string arg = "(?:" + name + ")|(?:" + flt + ")";
    // note: this includes the surrounding (...) that are needed to capture the values
    const std::string carg = ws + "((?:" + name + ")|(?:" + flt + "))" + ws;
    std::string args;
    for (unsigned i=0; i < arity; ++i)
    {
        if (i!=0)
        {
            args += ',';
        }
        args += carg;
    }
    const std::string cfname = '(' + name + ')';
    const std::string popen = "\\(";
    const std::string pclose = "\\)";
    const std::string func_call = cfname + popen + args + pclose;

    /* 2. Set up a regular expression object and parse the string. Issue an
     *    error if the string cannot be parsed as a call to a function with
     *    'arity' arguments.
     */
    const std::regex r{func_call};
    std::smatch m;
    if (!std::regex_match(str, m, r))
    {
        throw std::runtime_error("Error parsing '" + str + "' as a function with " + std::to_string(arity) + " arguments.");
    }
    // m contains the full string, followed by the function name and the arguments
    if (m.size()!=arity+2)
    {
        throw std::logic_error("INTERNAL ERROR: parseFunctionCall: unexpected number of arguments parsed.");
    }

    /* 3. Copy the function name and arguments into data members. The latter
     *    are converted to json objects. When an argument is not a valid JSON
     *    value, we assume it is a string (an identifier).
     */
    m_name = m[1];
    m_args.resize(arity);
    for (unsigned i=0; i!=arity; ++i)
    {
        const std::string& s = m[i+2];
        try
        {
            m_args[i] = json_type::parse(s.begin(),s.end());
        }
        catch(std::exception& exc)
        {
            // conversion to a JSON value failed. Assume it is a string.
            m_args[i] = s;
        }
    }
}

json_type makeJsonFromFunctionCall(
    const std::string& str,
    const std::string& fkey,
    const std::vector<std::string>& akeys)
{
    FunctionCall fcall(akeys.size(),str);
    json_type res;
    res[fkey] = fcall.name();
    for (unsigned i=0; i!= fcall.arity(); ++i)
    {
        res[akeys[i]] = fcall.arg(i);
    }
    return res;
}


bool getValue(const std::string &valueString, double &value)
{
    std::stringstream ss{valueString};
    return (ss >> value) && ss.eof();
}

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

std::string searchAndReplaceCopy(
                            const std::string& str,
                            const std::string& key,
                            const std::string& repl)
{
        std::string copy(str);
        return searchAndReplaceInPlace(copy,key,repl);
}

bool removeComments(std::istream& is, std::string &dest)
{
    dest.clear();
    std::string line;
    while (!is.eof())
    {
        std::getline(is,line);
        const std::string::size_type comment_pos = line.find('%');
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

bool stringBufferFromFile(const std::filesystem::path &fileName, std::string &dest)
{
    std::ifstream is(fileName);
    return is && removeComments(is,dest);
}

} // namespace Parse
} // namespace loki
