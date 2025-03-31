/** \file
 *
 *  Implementations of classes that represent terms of the power balance.
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

#include "LoKI-B/Power.h"

namespace loki
{

PowerTerm::PowerTerm()
: forward{0.0}, backward{0.0}
{
}

PowerTerm& PowerTerm::operator+=(const PowerTerm& term)
{
    forward += term.forward;
    backward += term.backward;
    return *this;
}

GasPower::GasPower()
{
}

GasPower& GasPower::operator+=(const GasPower& term)
{
    excitation += term.excitation;
    vibrational += term.vibrational;
    rotational += term.rotational;
    ionization += term.ionization;
    attachment += term.attachment;
    return *this;
}

Power::Power()
{
    clear();
}

void Power::clear()
{
    field = 0.0;
    elasticNet = 0.0;
    elasticGain = 0.0;
    elasticLoss = 0.0;
    carNet = 0.0;
    carGain = 0.0;
    carLoss = 0.0;

    inelastic = 0.0;
    superelastic = 0.0;
    eDensGrowth = 0.0;
    electronElectronNet = 0.0;
    electronElectronGain = 0.0;
    electronElectronLoss = 0.0;
    balance = 0.0;
    relativeBalance = 0.0;
    reference = 0.0;
}

Power& Power::operator+=(const GasPower &gasPower)
{
    // invoke += on our base class
    GasPower::operator+=(gasPower);
    return *this;
}

} // namespace loki
