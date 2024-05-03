/** \file
 *
 *  Code for reading gas properties and storing the results in a JSON object.
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
 *  \author Jan van Dijk
 *  \date   2. May 2024
 */

#ifndef LOKI_CPP_GASPROPERTIES_H
#define LOKI_CPP_GASPROPERTIES_H

#include "LoKI-B/json.h"
#include <filesystem>

namespace loki {

/** This class manages the information that is contained in a LoKI-B
 *  'gasProperties' section. It stores the data in a private json object.
 *  The elements' keys are the property names, the values are JSON objects.
 *  The latter's elements' keys are the names of the gases for which the
 *  property is available, and the value of that property.
 *
 *  \author Jan van Dijk
 *  \date   2. May 2024
 */
class GasProperties
{
public:
    GasProperties();
    GasProperties(const std::filesystem::path &basePath, const json_type& pnode);
    const json_type& data() const { return m_data; }
    json_type& data() { return m_data; }
private:
    json_type readGasPropertyFile(const std::filesystem::path& fname) const;
    json_type m_data;
};

} // namespace loki

#endif // LOKI_CPP_GASPROPERTIES_H
