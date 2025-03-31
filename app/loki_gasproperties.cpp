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
 *  \date   2 May 2024
 */

#include "LoKI-B/GasProperties.h"
#include "LoKI-B/LegacyToJSON.h"
#include "LoKI-B/json.h"
#include <stdexcept>

int main(int argc, const char *argv[])
try
{
    if (argc <3 || argc > 5)
    {
        throw std::runtime_error("Usage: loki_gasproperies <input file> <json pointer> [<property name> [<gas name>]].\n"
                                 "This program reads the gasProperties section that is defined by the input file name\n"
                                 "and json pointer. It prints the known properties to std::cout. If a third argument\n"
                                 "is provided, only the information about the requested property is listed. If a fourth\n"
                                 "argument is provided, the property will only be shown for te given gas.\n"
                                 "\n"
                                 "Example: ./app/loki_gasproperties ../input/default_lokib_setup.in /electronKinetics/gasProperties mass N4\n"
                                 "will print the mass of N4, when argument N4 is omitted, all known masses are printed.\n");
    }
    loki::json_type json;
    const std::filesystem::path path(argv[1]);
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
    const loki::json_type::json_pointer ptr{argv[2]};
    const loki::GasProperties gasProps(path,json[ptr]);
    switch (argc)
    {
        case 3:
            std::cout << gasProps.data().dump(2) << std::endl;
        break;
        case 4:
            std::cout << gasProps.data().at(argv[3]).dump(2) << std::endl;
        break;
        case 5:
            std::cout << gasProps.data().at(argv[3]).at(argv[4]).dump(2) << std::endl;
        break;
    }
    return 0;
}
catch (std::exception &exc)
{
    std::cerr << exc.what() << std::endl;
    return 1;
}
