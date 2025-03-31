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

#include "LoKI-B/json.h"
#include "LoKI-B/Log.h"
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace loki
{

json_type read_json_from_stream(std::istream &is)
{
    return json_type::parse(is);
}

json_type read_json_from_file(const std::filesystem::path &fname)
{
    std::ifstream is(fname);
    if (!is)
    {
        Log<Message>::Error("Could not open file '", fname.generic_string(), "' for reading.");
    }
    try
    {
        return read_json_from_stream(is);
    }
    catch (std::exception &exc)
    {
        throw std::runtime_error("Error reading file '" + fname.generic_string() + "':\n" + std::string(exc.what()));
    }
}

} // namespace loki
