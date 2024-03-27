/** \file
 *
 *  Implementation of a lookup table class for LoKI-B
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
 *  \author Jan van Dijk and Daan Boer
 *  \date November 2020
 */

#include "LoKI-B/LookupTable.h"
#include <stdexcept>

namespace loki {

LookupTable::LookupTable(const Vector& x, const Vector& y)
    : m_x(x), m_y(y)
{
    if (m_x.size()!=m_y.size())
    {
        throw std::runtime_error("Dimensions of abscissa and ordinate do not match.");
    }
    if (m_x.size()<2)
    {
        throw std::runtime_error("At least two points are required.");
    }
    for (Index i=1; i!=m_x.size(); ++i)
    {
        if (m_x[i]<=m_x[i-1])
        {
            throw std::runtime_error("Abscissa values must appear in ascending order.");
        }
    }
}

LookupTable LookupTable::create(const json_type& data)
{
    using PairVector = std::vector<std::pair<double, double>>;
    const PairVector pairs(data.at("values").get<PairVector>());
    Vector x,y;
    x.resize(pairs.size());
    y.resize(pairs.size());
    for (PairVector::size_type i = 0; i != pairs.size(); ++i)
    {
        x[i] = pairs[i].first;
        y[i] = pairs[i].second;
    }
    return LookupTable(x,y);
}

LookupTable LookupTable::create(std::istream& is)
{
    std::vector<double> x, y;
    std::string line;

    while (std::getline(is, line))
    {
        if (line.substr(0, 2) == "--")
        {
            break;
        }
    }
    while (std::getline(is, line))
    {
        if (line.size() && line[0]=='#')
        {
            continue;
        }
        std::stringstream ss(line);
        double energy, value;
        ss >> energy >> value;
        if (!ss)
        {
            // This line did not start with two valid numbers.
            // Check that we indeed reached the end-of-table
            // marker "--" (there may be more "-").
            if (line.substr(0, 2) != "--")
            {
                throw std::runtime_error("Bad data line '" + line + "': expected numbers or '--'");
            }
            // OK, -- was found. Stop reading.
            break;
        }
        x.emplace_back(energy);
        y.emplace_back(value);
    }
    return LookupTable(
            Vector::Map(x.data(), x.size()),
            Vector::Map(y.data(), y.size())
    );
}

} // namespace loki
