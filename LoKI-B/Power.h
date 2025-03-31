/** \file
 *
 *  Declaration of classes that represent terms of the power balance.
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
    PowerTerm();
    double net() const { return forward + backward; }
    PowerTerm& operator+=(const PowerTerm& term);

    double forward;
    double backward;
};

/** While ionization, attachment only have forward 'inE' parts, we represent
 *  these terms with a full PowerTerm object. That hardly costs anything, and
 *  has the advantage of uniformity and nicer semantics: if somebody is
 *  interested in the net ionization energy losses, he can now say
 *  ionization.net(). Of course this will still be equal to ionization.forward
 *  in the absence of reverse processes.
 */
struct GasPower
{
    GasPower();
    GasPower& operator+=(const GasPower& term);

    PowerTerm excitation;
    PowerTerm vibrational;
    PowerTerm rotational;
    PowerTerm ionization;
    PowerTerm attachment;
};

struct Power : GasPower
{
    Power();
    void clear();
    Power& operator+=(const GasPower &gasPower);

    double field;
    double elasticNet;
    double elasticGain;
    double elasticLoss;
    double carNet;
    double carGain;
    double carLoss;
    double inelastic;
    double superelastic;
    double eDensGrowth;
    /** Electron-electron collisions result in a redistribution of energy
     *  among the electrons, so this term should be zero. The value of this
     *  term is used as a convergence test in the code.
     */
    double electronElectronNet;
    double electronElectronGain;
    double electronElectronLoss;
    double balance;
    double relativeBalance;
    double reference;
};

} // namespace loki

#endif // LOKI_CPP_POWER_H
