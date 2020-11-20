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

bool stringBufferFromFile(const std::string &fileName, std::string &dest)
{
    std::ifstream is(INPUT "/" + fileName);
    return is && removeComments(is,dest);
}

} // namespace Parse
} // namespace loki
