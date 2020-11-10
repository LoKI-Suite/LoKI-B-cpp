//
// Created by daan on 15-5-19.
//

#include "LoKI-B/CrossSection.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/Parse.h"
#include <utility>

namespace loki
{

CrossSection::CrossSection(Grid *energyGrid, bool isElasticOrEffective, const json_type &cnf)
    : m_threshold(cnf.contains("threshold") ? cnf.at("threshold").get<double>() : 0.0),
    m_energyGrid(energyGrid),
    m_isElasticOrEffective(isElasticOrEffective),
    m_lut(LookupTable::create(cnf))
{
    this->interpolate();
    energyGrid->updatedMaxEnergy1.addListener(&CrossSection::interpolate, this);
}

CrossSection::CrossSection(double threshold, Grid *energyGrid, bool isElasticOrEffective, std::istream &in)
    : m_threshold(threshold),
    m_energyGrid(energyGrid),
    m_isElasticOrEffective(isElasticOrEffective),
    m_lut(LookupTable::create(in))
{
    this->interpolate();
    energyGrid->updatedMaxEnergy1.addListener(&CrossSection::interpolate, this);
}

CrossSection::CrossSection(double threshold, Grid *energyGrid, bool isElasticOrEffective, Vector rawEnergyData,
                           Vector rawCrossSection)
    : m_threshold(threshold),
    m_energyGrid(energyGrid),
    m_isElasticOrEffective(isElasticOrEffective),
    m_lut(rawEnergyData,rawCrossSection)
{
    this->interpolate();
    energyGrid->updatedMaxEnergy1.addListener(&CrossSection::interpolate, this);
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
    /** \todo Could we also use if(threshold<=0) instead of if (!m_isElasticOrEffective)?
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
        /* 0.0, 0.0: these are the values that are assumed if the enery is not
         * in the range of the table.
         */
	result[gridIndex] = m_lut.interpolate_or_set(energies[gridIndex],0.0,0.0);
    }
}

} // namespace loki
