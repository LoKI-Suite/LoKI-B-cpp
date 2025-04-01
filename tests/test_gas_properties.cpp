/** \file
 *
 *  Tests for function offSideToJSON.
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
 *  \date   3 May 2024
 */

#include "LoKI-B/GasProperties.h"
#include "LoKI-B/Constant.h"
#include "tests/TestUtilities.h"
#include <iostream>
#include <sstream>

void test_manual()
try
{
    loki::GasProperties gasProps;

    // gasProps does not have property "rotationalConstant" by default
    test_expr(gasProps.has("rotationalConstant")==false);

    // the electron mass is added by default
    test_expr(gasProps.has("mass")==true);
    test_expr(gasProps.has("mass","e")==true);
    test_expr(gasProps.get("mass","e")==loki::Constant::electronMass);
    test_expr(gasProps.get("mass","e",1)==loki::Constant::electronMass);

    test_expr(gasProps.has("mass","X")==false);
    try
    {
        test_expr(gasProps.get("mass","X")==loki::Constant::electronMass);
        std::cout << "ERROR: expected that an exception would be thrown." << std::endl;
        ++nerrors;
    }
    catch (std::exception& exc)
    {
        // OK. An exception was expected
    }
    // test that the alternative value is returned when the property is not available
    test_expr(gasProps.get("mass","X",42)==42);

    // add the rotationalConstant of gas "X"
    // NOTE: non-sensical numbers are used here for testing purposes.
    gasProps.set("rotationalConstant","X",1);
    test_expr(gasProps.has("rotationalConstant")==true);
    test_expr(gasProps.has("rotationalConstant","X")==true);
    test_expr(gasProps.get("rotationalConstant","X")==1);
}
catch(std::exception& exc)
{
    std::cout << "ERROR: unexpected exception: " << exc.what() << std::endl;
    ++nerrors;
}

int main()
{
    test_manual();

    test_report;
    return nerrors;
}
