/** \file
 *
 *  Structures that hold the swarm parameters and rate coefficients.
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
 *  \author Daan Boer and
 *  \date   2 July 2019
 */

#ifndef LOKI_CPP_MACROSCOPICQUANTITIES_H
#define LOKI_CPP_MACROSCOPICQUANTITIES_H

#include <complex>

namespace loki {

struct SwarmParameters
{
    double redDiffCoeff{0.};
    double redMobCoeff{0.};
    double redTownsendCoeff{0.};
    double redAttCoeff{0.};
    double redDiffCoeffEnergy{0.};
    double redMobilityEnergy{0.};
    std::complex<double> redMobilityHF{0.,0.};
    double meanEnergy{0.};
    double characEnergy{0.};
    double Te{0.};
    double driftVelocity{0.};
};

// Forward declaration
class EedfCollision;

struct RateCoefficient
{
    const EedfCollision *collision;
    double inelastic;
    double superelastic;
};

} // namespace loki

#endif // LOKI_CPP_MACROSCOPICQUANTITIES_H
