/** \file
 *
 *  A demonstration of function fgPrimeEnergyIntegral of the LoKI-B project.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2022 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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
 *  Test functions energyIntegral and fgPrimeEnergyIntegral.
 *  energyIntegral approximates \f$ \int_0^{u_max} f(u)g(u) du \f$;
 *  fgPrimeEnergyIntegral approximates \f$ \int_0^{u_max} f*(dg/du) du.
 *  We test these function for some stylized functions f and g for which
 *  the result is known.
 *
 *  \author Jan van Dijk
 *  \date   7 June 2024
 */

#include "LoKI-B/Grid.h"
#include "LoKI-B/Integrals.h"
#include "LoKI-B/Operators.h"
#include "tests/TestUtilities.h"

using namespace loki;

/* f = 1, g=u =>
 * int_0^{u_max} f(u)(dg/du) du = u_max.
 * int_0^{u_max} f(u)g(u) du    = u_max^2/2.
 */
void test1()
{
    const Grid grid(10,10);
    Vector f(grid.nCells()); f.fill(1);
    const Vector g = grid.getCells();
    test_expr(std::abs(energyIntegral(grid,f,g) - grid.uMax()*grid.uMax()/2)
        < grid.nCells()*std::numeric_limits<double>::epsilon());
    test_expr(std::abs(fgPrimeEnergyIntegral(grid,f,g) - grid.uMax())
        < grid.nCells()*std::numeric_limits<double>::epsilon());
}

/* f = u, g=u =>
 * int_0^{u_max} f(u)(dg/du) du = u_max^2/2.
 * int_0^{u_max} f(u)g(u) du    = u_max^3/3.
 */
void test2()
{
    const Grid grid(10,10);
    const Vector f = grid.getCells();
    const Vector g = grid.getCells();
    /** \todo Establish an upper bound for the error, given n and the
     *  discretisation scheme, and use that:
     */
    test_expr(std::abs(energyIntegral(grid,f,g) - grid.uMax()*grid.uMax()*grid.uMax()/3)
        < 1);
    test_expr(std::abs(fgPrimeEnergyIntegral(grid,f,g) - grid.uMax()*grid.uMax()/2)
        < 10*grid.nCells()*std::numeric_limits<double>::epsilon());
}

int main()
{
    test1();
    test2();
    test_report; 
    return nerrors; 
}

