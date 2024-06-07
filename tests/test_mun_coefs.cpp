/** \file
 *
 *  A demonstration of function calculateMuNCoefs of the LoKI-B project.
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
 *  Test the calculation of the coefficients that are used to calculate
 *  the electron mobility for a given EEDF.
 *
 *  \author Jan van Dijk
 *  \date   7 June 2024
 */

#include "LoKI-B/Constant.h"
#include "LoKI-B/Grid.h"
#include "LoKI-B/Operators.h"
#include "tests/TestUtilities.h"

using namespace loki;

/* The coefficients muNCoefs that are returned by calculateMuNCoefs are such
 * that muNCoefs.dot(f) \approx gamma*int_0^{u_max} D0*(df/du) du. We test
 * these coefficients for some stylized functions D0 and f for which the
 * integration should be exact.
 */

// D0 = 1/gamma, f=u => int_0^{u_max} (df/du)du = u_max.
void test1()
{
    const Grid grid(10,10);
    Vector D0(grid.nCells()); D0.fill(1/SI::gamma);
    const Vector eedf = grid.getCells();
    const Vector muNCoefs(calculateMuNCoefs(D0));
    test_expr(std::abs(muNCoefs.dot(eedf)-grid.uMax())
        < grid.nCells()*std::numeric_limits<double>::epsilon());
}

//  D0 = u/gamma, f=u => int_0^{u_max} u(df/du)du = u_max^2/2.
void test2()
{
    const Grid grid(10,10);
    Vector D0 = grid.getCells()/SI::gamma;
    const Vector eedf = grid.getCells();
    const Vector muNCoefs(calculateMuNCoefs(D0));
    test_expr(std::abs(muNCoefs.dot(eedf)-grid.uMax()*grid.uMax()/2)
        < grid.nCells()*std::numeric_limits<double>::epsilon());
}

int main()
{
    test1();
    test2();
    test_report; 
    return nerrors; 
}

