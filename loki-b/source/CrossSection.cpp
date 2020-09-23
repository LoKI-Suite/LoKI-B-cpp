#include <utility>

//
// Created by daan on 15-5-19.
//

#include "LoKI-B/CrossSection.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/Parse.h"

namespace loki
{

CrossSection::CrossSection(Grid *energyGrid, bool isElasticOrEffective, const json_type &cnf)
    : threshold(cnf.contains("threshold") ? cnf.at("threshold").get<double>() : 0.0), energyGrid(energyGrid),
      isElasticOrEffective(isElasticOrEffective)
{

    using PairVector = std::vector<std::pair<double, double>>;
    const PairVector tmp(cnf.at("data"));
    rawEnergyData.resize(tmp.size());
    rawCrossSection.resize(tmp.size());
    for (PairVector::size_type i = 0; i != tmp.size(); ++i)
    {
        rawEnergyData[i] = tmp[i].first;
        rawCrossSection[i] = tmp[i].second;
    }
#if 0
std::cout << "isElasticOrEffective: " << isElasticOrEffective << std::endl;
std::cout << "THRESHOLD: " << threshold << std::endl;
for (unsigned i=0; i!=rawEnergyData.size(); ++i)
{
  std::cout << rawEnergyData[i] << " " << rawCrossSection[i] << std::endl;
}
#endif
    this->interpolate();
    this->energyGrid->updatedMaxEnergy1.addListener(&CrossSection::interpolate, this);
}

//    CrossSection::CrossSection(const double threshold, Grid *energyGrid, bool isElasticOrEffective)
//            : threshold(threshold), energyGrid(energyGrid), isElasticOrEffective(isElasticOrEffective) {}

CrossSection::CrossSection(double threshold, Grid *energyGrid, bool isElasticOrEffective, std::istream &in)
    : threshold(threshold), energyGrid(energyGrid), isElasticOrEffective(isElasticOrEffective)
{
    std::vector<double> rawEnergyVector, rawCrossSectionVector;

    Parse::rawCrossSectionFromStream(rawEnergyVector, rawCrossSectionVector, in);

    rawEnergyData = Vector::Map(rawEnergyVector.data(), rawEnergyVector.size());
    rawCrossSection = Vector::Map(rawCrossSectionVector.data(), rawCrossSectionVector.size());
#if 0
std::cout << "isElasticOrEffective: " << isElasticOrEffective << std::endl;
std::cout << "THRESHOLD: " << threshold << std::endl;
for (unsigned i=0; i!=rawEnergyData.size(); ++i)
{
  std::cout << rawEnergyData[i] << " " << rawCrossSection[i] << std::endl;
}
#endif

    this->interpolate();
    this->energyGrid->updatedMaxEnergy1.addListener(&CrossSection::interpolate, this);
}

CrossSection::CrossSection(double threshold, Grid *energyGrid, bool isElasticOrEffective, Vector rawEnergyData,
                           Vector rawCrossSection)
    : threshold(threshold), energyGrid(energyGrid), isElasticOrEffective(isElasticOrEffective),
      rawEnergyData(std::move(rawEnergyData)), rawCrossSection(std::move(rawCrossSection))
{
    this->interpolate();
    this->energyGrid->updatedMaxEnergy1.addListener(&CrossSection::interpolate, this);
}

void CrossSection::interpolate()
{

    interpolate(energyGrid->getNodes(), *this);
}

void CrossSection::interpolate(const Vector &energies, Vector &result)
{
    const auto &gridSize = energies.size();

    result.resize(gridSize);
    result.setZero(gridSize);

    uint32_t csIndex = 0, gridIndex = 0;

    if (!isElasticOrEffective)
    {
        for (uint32_t i = 0; i < gridSize; ++i)
        {
            if (energies[i] > threshold)
            {
                gridIndex = i;
                break;
            }
        }
    }

    for (; gridIndex < gridSize; ++gridIndex)
    {
        while (csIndex < rawCrossSection.size() && rawEnergyData[csIndex] < energies[gridIndex])
        {

            ++csIndex;
        }

        if (csIndex >= rawCrossSection.size())
        {
            //                (*this)[gridIndex] = rawCrossSection.back().second;
            result[gridIndex] = 0.;
            continue;
        }

        if (csIndex == 0)
        {
            if (rawEnergyData[csIndex] == energies[gridIndex])
            {
                result[gridIndex] = rawCrossSection[csIndex];
            }
            else
            {
                result[gridIndex] = 0.;
            }
            continue;
        }

        const double prevEnergy = rawEnergyData[csIndex - 1];
        const double nextEnergy = rawEnergyData[csIndex];
        const double prevCS = rawCrossSection[csIndex - 1];
        const double nextCS = rawCrossSection[csIndex];
        const double alpha = (energies[gridIndex] - prevEnergy) / (nextEnergy - prevEnergy);

        result[gridIndex] = (1. - alpha) * prevCS + alpha * nextCS;
    }
}

Vector &CrossSection::raw()
{
    return rawCrossSection;
}

Vector &CrossSection::energies()
{
    return rawEnergyData;
}

const Grid *CrossSection::getGrid() const
{
    return energyGrid;
}

} // namespace loki
