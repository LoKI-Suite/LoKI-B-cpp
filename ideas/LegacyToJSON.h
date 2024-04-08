/** \file
 *
 *  Convert legacy LoKI-B input files to JSON.
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
 *  \date   21. May 2019
 */

#ifndef LOKI_CPP_LEGACYTOJSON_H
#define LOKI_CPP_LEGACYTOJSON_H

#include "LoKI-B/json.h"
#include <iostream>
#include <filesystem>

namespace loki
{

json_type legacyToJSON(std::istream& is);
json_type legacyToJSON(const std::filesystem::path& fname);

} // namespace loki



#endif // LOKI_CPP_LEGACYTOJSON_H

