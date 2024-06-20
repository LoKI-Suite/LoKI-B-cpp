/** \file
 *
 *  Declaration of classes that represent terms of the power balance.
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
 *  \author Daan Boer
 *  \date   24 June 2019
 */

#ifndef LOKI_CPP_POWER_H
#define LOKI_CPP_POWER_H

#include "LoKI-B/Enumeration.h"

namespace loki
{

struct PowerTerm
{
    PowerTerm() : forward{0.0}, backward{0.0} { }
    double net() const { return forward + backward; }
    PowerTerm& operator+=(const PowerTerm& term)
    {
        forward += term.forward;
        backward += term.backward;
        return *this;
    }
    double forward;
    double backward;
};

/** \todo Document that this was changed during the development of the C++ code:
 *  while ionization, attachment only have forward 'inE' parts, we represent these
 *  terms with a full PowerTerm object. That hardly costs anything, and has the
 *  advantage of uniformity and nicer semantics: if somebody is interested in the
 *  net ionization energy losses, he can now say ionization.net(). Of course this
 *  will still be equal to ionization.forward in the absenace of reverse processes.
 */
struct GasPower
{
    GasPower() { }
    PowerTerm excitation;
    PowerTerm vibrational;
    PowerTerm rotational;
    PowerTerm ionization;
    PowerTerm attachment;
    GasPower& operator+=(const GasPower& term)
    {
        excitation += term.excitation;
        vibrational += term.vibrational;
        rotational += term.rotational;
        ionization += term.ionization;
        attachment += term.attachment;
        return *this;
    }
};

struct Power : GasPower
{
    double field{0.};
    double elasticNet{0.}, elasticGain{0.}, elasticLoss{0.};
    double carNet{0.}, carGain{0.}, carLoss{0.};

    double inelastic{0.};
    double superelastic{0.};
    double eDensGrowth{0.};
    /** \todo What is the meaning of this term? This is a redistribution of
     *        energy among the electrons, so I would think that this should
     *        be just zero. Is this meant as a check?
     */
    double electronElectronNet{0.};
    double electronElectronGain{0.};
    double electronElectronLoss{0.};
    double balance{0.};
    double relativeBalance{0.};
    double reference{0.};

    Power& operator+=(const GasPower &gasPower)
    {
	// invoke += on our base class
        GasPower::operator+=(gasPower);
	return *this;
    }
};

} // namespace loki

#endif // LOKI_CPP_POWER_H
