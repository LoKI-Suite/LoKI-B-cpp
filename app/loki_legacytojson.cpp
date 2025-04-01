/** \file A utility for converting a legacy LoKI-B file to LoKI-B JSON format.
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
 *  \author Jan van Dijk
 *  \date   30 April 2024
 */

#include "LoKI-B/LegacyToJSON.h"
#include "LoKI-B/json.h"
#include <stdexcept>

int main(int argc, const char *argv[])
try
{
    loki::json_type json;
    switch (argc)
    {
    case 1:
        json = loki::legacyToJSON(std::cin);
        break;
    case 2: {
        const auto path = std::filesystem::path(argv[1]);
        if (path.extension() == ".json")
        {
            json = loki::legacyToJSON(loki::read_json_from_file(path));
        }
        else if (path.extension() == ".in")
        {
            json = loki::legacyToJSON(std::filesystem::path(argv[1]));
        }
        else
        {
            throw ::std::runtime_error("Expected a file with extension '.json', or '.in'.");
        }
        break;
    }
    default:
        throw std::runtime_error("Usage: loki_legacytojson [input file]\n"
                                 "This program converts a legacy LoKI-B file to json format\n"
                                 "and makes the necessary modifications to allow the usage of\n"
                                 "this file with LoKI-B. The resulting JSON object is printed\n"
                                 "to the console. If no input file is provided, input will be\n"
                                 "read from the standard input stream.");
        return 1;
        break;
    }
    std::cout << json.dump(2) << std::endl;
    return 0;
}
catch (std::exception &exc)
{
    std::cerr << exc.what() << std::endl;
    return 1;
}
