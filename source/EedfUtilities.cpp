/** \file
 *
 *  Implementation of LoKI-B EEDF-related utility functions.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2024 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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
 *  \date   29 July 2023
 */

#include "LoKI-B/EedfUtilities.h"
#include "LoKI-B/GridOps.h"
#include <cmath>

namespace loki {

double getMeanEnergy(const Vector& edf, const Grid& grid)
{
    return energyIntegral(grid,grid.getCells().cwiseProduct(grid.getCells().cwiseSqrt()),edf);
}

void normalizeEDF(Vector& edf, const Grid& grid)
{
    edf /= energyIntegral(grid,grid.getCells().cwiseSqrt(),edf);
}

void makePrescribedEDF(Vector& edf, const Grid& grid, double s, double T_eV, bool normalize)
{
    edf.resize(grid.nCells());
    const double gamma_3_2s = std::tgamma(3/(2*s));
    const double gamma_5_2s = std::tgamma(5/(2*s));
    const double factor = s*std::pow(gamma_5_2s,1.5) * std::pow(gamma_3_2s,-2.5) * std::pow(2/(3*T_eV),1.5);
    for (Grid::Index k = 0; k != grid.nCells(); ++k)
    {
        const double energy = grid.getCell(k);
        edf(k) = factor * std::exp(-std::pow(energy*gamma_5_2s*2/(gamma_3_2s*3*T_eV),s));
    }
    if (normalize)
    {
        normalizeEDF(edf,grid);
    }
}

Vector makePrescribedEDF(const Grid& grid, double s, double T_eV, bool normalize)
{
    Vector edf(grid.nCells());
    makePrescribedEDF(edf,grid,s,T_eV,normalize);
    return edf;
}

} // namespace loki

