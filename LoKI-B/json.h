/** \file
 *
 *  Support for nlohmann_json.
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
 *  \author Jan van Dijk (C++ version)
 *  \date   July 2020
 */

#ifndef LOKI_CPP_JSON_H
#define LOKI_CPP_JSON_H

/** Instruct nlohmann::json objects to keep track of their parents, so
 *  sensible error messages can be produced when (for example) access
 *  to an element or a type conversion fails.
 *
 *  See https://json.nlohmann.me/api/macros/json_diagnostics/
 *
 *  gcc is known to produce false positives when the -Warray-bounds
 *  flag is used --- and indeed it does. The warning messages are so
 *  noisy that the compiler output is rendered useless. That flag has
 *  therefore been disabled in CMakeLists.txt. For more details, see
 *  https://github.com/nlohmann/json/issues/3808
 */
#define JSON_DIAGNOSTICS 1

#include <nlohmann/json.hpp>
#include "LoKI-B/Exports.h"
#include <filesystem>
#include <iostream>
#include <string>

namespace loki
{

using json_type = nlohmann::json;

lokib_export json_type read_json_from_stream(std::istream &is);
lokib_export json_type read_json_from_file(const std::filesystem::path &fname);

} // namespace loki

#endif // LOKI_CPP_JSON_H
