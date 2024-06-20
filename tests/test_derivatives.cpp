/** \file
 *
 *  A demonstration of function cellDerivative of the LoKI-B project.
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
 *  Test function cellDerivative, which returns an approximation of a
 *  field that is defined in the cells of a grid in (again) the cells.
 *
 *  \author Jan van Dijk
 *  \date   14 June 2024
 */

#include "LoKI-B/Grid.h"
#include "LoKI-B/GridOps.h"
#include "tests/TestUtilities.h"

using namespace loki;

/* Test that u_face is interpolated to u_cell.
 */
void test1(const Grid& grid)
{
    const Vector uprime = cellDerivative(grid,grid.getCells());
    Vector ones(grid.nCells()); ones.fill(1.0);
    test_expr ( (uprime-ones).norm() < grid.nCells()*std::numeric_limits<double>::epsilon() );
}

int main()
{
    const Grid grid1(10,10);
    Vector nodes(5); nodes << 0.0, 0.5, 0.75, 0.875, 1.0;
    const Grid grid2(nodes,10);

    test1(grid1);
    test1(grid2);
    test_report; 
    return nerrors; 
}

