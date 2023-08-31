/** \file
 *
 *  A demonstration of class ElectronElectronOperator of the LoKI-B project.
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
 *  Test the calculation of the (reduced) particle and energy mobility
 *  and diffusion coefficients of the electrons.
 *
 *  \author Jan van Dijk
 *  \date   August 2023
 */

#include "LoKI-B/Grid.h"
#include "LoKI-B/EedfUtilities.h"
#include "LoKI-B/Constant.h"
#include <iostream>
#include <fstream>

#include <Eigen/SparseLU>
#include <Eigen/OrderingMethods>

using namespace loki;

double calculateDN(const Grid& grid, const Vector& D0, const Vector& eedf)
{
    // This is 33a from \cite Manual_1_0_0
    return SI::gamma * grid.du() * D0.dot(eedf);
}

double calculateMuNOld(const Grid& grid, const Vector& D0, const Vector& eedf)
{
    Vector U0sup(grid.nCells());
    Vector U0inf(grid.nCells());
    U0sup[0] = 0.;
    U0inf[grid.nCells() - 1] = 0.;
    for (Grid::Index j = 0; j < grid.nCells(); ++j)
    {
        if (j != 0)
            U0sup[j] = +1./(2. * grid.du()) * D0[j - 1];

        if (j != grid.nCells() - 1)
            U0inf[j] = -1./(2. * grid.du()) * D0[j + 1];
    }
    const Vector U0 = U0sup + U0inf;

    /* This is 33b from \cite Manual_1_0_0, multiplied with -1.
     */
    return - SI::gamma * grid.du() * U0.dot(eedf);
}

double calculateMuNNew(const Grid& grid, const Vector& D0, const Vector& eedf)
{
    Vector dDdu_l(grid.getCells().size());
    Vector dDdu_h(grid.getCells().size());
    dDdu_l[0] = -D0[0]/grid.du();
    dDdu_h[0] = +D0[1]/grid.du();
    for (Grid::Index c=1; c!=grid.getCells().size()-1; ++c)
    {
        dDdu_l[c] = - D0[c-1]/(2*grid.du());
        dDdu_h[c] = + D0[c+1]/(2*grid.du());
    }
    const Grid::Index nc = grid.getCells().size();
    dDdu_l[nc-1] = -D0[nc-2]/grid.du();
    dDdu_h[nc-1] = +D0[nc-1]/grid.du();
    return SI::gamma * grid.du() * (dDdu_l + dDdu_h).dot(eedf);
}

int main(int argc, const char* argv[])
{
    const unsigned nCells = 1000;
    const double uMax = 10; // eV

    const double shape = 1.0; // 1: Maxwell, 2: Druyvesteyn

    const Grid grid(nCells,uMax);
    Vector cellTotalCrossSection(grid.getCells().size());
    for (Grid::Index c=0; c!=grid.getCells().size(); ++c)
    {
        // in m^2
        cellTotalCrossSection[c] = 1e-20;
        //cellTotalCrossSection[c] = 1/std::sqrt(grid.getCells()[c]);
    }
    std::cout << "# TeV\tDN\tmuN(old)\tEinstein(old)\tmuN(new)\tEinstein(new)" << std::endl;
    for (double Te =300; Te<=10000; Te+=100)
    {
        const double TeV = Te*Constant::boltzmann/Constant::electronCharge;
        const Vector eedf = makePrescribedEDF(grid,shape,TeV);
        const Vector D0 = grid.getCells().array() / (3. * cellTotalCrossSection.array());
        const double DN = calculateDN(grid,D0,eedf);
        const double muNOld = calculateMuNOld(grid,D0,eedf);
        const double muNNew = calculateMuNNew(grid,D0,eedf);
        std::cout << TeV
            << '\t' << DN
            << '\t' << muNOld
            << '\t' << Constant::boltzmann*Te*muNOld/(DN*Constant::electronCharge)
            << '\t' << muNNew
            << '\t' << Constant::boltzmann*Te*muNNew/(DN*Constant::electronCharge)
            << std::endl;
            // note: columns 4 and 6 should be 1 for Maxwell
    }
    return 0;
}
