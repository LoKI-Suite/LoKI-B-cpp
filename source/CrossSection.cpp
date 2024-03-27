/** \file
 *
 *  Implementation of a class that represents cross sections.
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
 *  \author Daan Boer and Jan van Dijk (C++ version)
 *  \date   15 May 2019
 */

#include "LoKI-B/CrossSection.h"
#include "LoKI-B/Log.h"
#include <utility>

namespace loki
{

CrossSection::CrossSection(const Grid *energyGrid, bool isElasticOrEffective, const json_type &cnf)
    : m_threshold(cnf.contains("threshold") ? cnf.at("threshold").get<double>() : 0.0),
    m_energyGrid(energyGrid),
    m_isElasticOrEffective(isElasticOrEffective),
    m_lut(LookupTable::create(cnf["data"]))
{
    this->interpolate();
    energyGrid->updatedMaxEnergy.addListener(&CrossSection::interpolate, this);
}

CrossSection::CrossSection(double threshold, const Grid *energyGrid, bool isElasticOrEffective, std::istream &in)
    : m_threshold(threshold),
    m_energyGrid(energyGrid),
    m_isElasticOrEffective(isElasticOrEffective),
    m_lut(LookupTable::create(in))
{
    this->interpolate();
    energyGrid->updatedMaxEnergy.addListener(&CrossSection::interpolate, this);
}

CrossSection::CrossSection(double threshold, const Grid *energyGrid, bool isElasticOrEffective, Vector rawEnergyData,
                           Vector rawCrossSection)
    : m_threshold(threshold),
    m_energyGrid(energyGrid),
    m_isElasticOrEffective(isElasticOrEffective),
    m_lut(rawEnergyData,rawCrossSection)
{
    this->interpolate();
    energyGrid->updatedMaxEnergy.addListener(&CrossSection::interpolate, this);
}

void CrossSection::interpolate()
{
    interpolate(m_energyGrid->getNodes(), *this);
}

void CrossSection::interpolate(const Vector &energies, Vector &result) const
{
    const Index nEnergies = energies.size();
    result.resize(nEnergies);
    Index gridIndex = 0;

    /** \todo Like in the old code, sigma is set to zero if the energy is
     *        exactly equal to the threshold. This choice should be documented.
     */
    /** \todo Could we also use if(threshold>0) instead of if (!m_isElasticOrEffective)?
     *        In that case, data member m_isElasticOrEffective and the constructor
     *        argument isElasticOrEffective) are no longer needed.
     */
    // if this is an inelastic process (with a threshold energy), set
    // the cross section values for energies up to the threshold to zero.
    if (!m_isElasticOrEffective)
    {
        for (; gridIndex!=nEnergies && energies[gridIndex]<=m_threshold; ++gridIndex)
        {
            result[gridIndex] = 0.0;
        }
    }
    for (; gridIndex != nEnergies; ++gridIndex)
    {
        /* 0.0, 0.0: these are the values that are assumed if the enery is
         * below or above the range of the table, respectively.
         */
	result[gridIndex] = m_lut.interpolate_or_set(energies[gridIndex],0.0,0.0);
    }
}

} // namespace loki
