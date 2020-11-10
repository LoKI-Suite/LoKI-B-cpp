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
    m_isElasticOrEffective(isElasticOrEffective)
{
    using PairVector = std::vector<std::pair<double, double>>;
    const PairVector tmp(cnf.at("data"));
    m_rawEnergyData.resize(tmp.size());
    m_rawCrossSection.resize(tmp.size());
    for (PairVector::size_type i = 0; i != tmp.size(); ++i)
    {
        m_rawEnergyData[i] = tmp[i].first;
        m_rawCrossSection[i] = tmp[i].second;
    }
    this->interpolate();
    energyGrid->updatedMaxEnergy1.addListener(&CrossSection::interpolate, this);
}

CrossSection::CrossSection(double threshold, Grid *energyGrid, bool isElasticOrEffective, std::istream &in)
    : m_threshold(threshold),
    m_energyGrid(energyGrid),
    m_isElasticOrEffective(isElasticOrEffective)
{
    std::vector<double> rawEnergyVector, rawCrossSectionVector;

    Parse::rawCrossSectionFromStream(rawEnergyVector, rawCrossSectionVector, in);

    m_rawEnergyData = Vector::Map(rawEnergyVector.data(), rawEnergyVector.size());
    m_rawCrossSection = Vector::Map(rawCrossSectionVector.data(), rawCrossSectionVector.size());

    this->interpolate();
    energyGrid->updatedMaxEnergy1.addListener(&CrossSection::interpolate, this);
}

CrossSection::CrossSection(double threshold, Grid *energyGrid, bool isElasticOrEffective, Vector rawEnergyData,
                           Vector rawCrossSection)
    : m_threshold(threshold),
    m_energyGrid(energyGrid),
    m_isElasticOrEffective(isElasticOrEffective),
    m_rawEnergyData(rawEnergyData),
    m_rawCrossSection(rawCrossSection)
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
    const Index gridSize = energies.size();

    result.resize(gridSize);
    result.setZero(gridSize);

    Index csIndex = 0, gridIndex = 0;

    if (!m_isElasticOrEffective)
    {
        for (Index i = 0; i < gridSize; ++i)
        {
            if (energies[i] > m_threshold)
            {
                gridIndex = i;
                break;
            }
        }
    }

    for (; gridIndex < gridSize; ++gridIndex)
    {
        while (csIndex < m_rawCrossSection.size() && m_rawEnergyData[csIndex] < energies[gridIndex])
        {
            ++csIndex;
        }

        if (csIndex >= m_rawCrossSection.size())
        {
            result[gridIndex] = 0.;
            continue;
        }

        if (csIndex == 0)
        {
            if (m_rawEnergyData[csIndex] == energies[gridIndex])
            {
                result[gridIndex] = m_rawCrossSection[csIndex];
            }
            else
            {
                result[gridIndex] = 0.;
            }
            continue;
        }

        const double prevEnergy = m_rawEnergyData[csIndex - 1];
        const double nextEnergy = m_rawEnergyData[csIndex];
        const double prevCS = m_rawCrossSection[csIndex - 1];
        const double nextCS = m_rawCrossSection[csIndex];
        const double alpha = (energies[gridIndex] - prevEnergy) / (nextEnergy - prevEnergy);

        result[gridIndex] = (1. - alpha) * prevCS + alpha * nextCS;
    }
}

} // namespace loki
