#include <utility>

//
// Created by daan on 15-5-19.
//

#include "CrossSection.h"
#include "Parse.h"
#include "Log.h"

namespace loki {

    CrossSection::CrossSection(const double threshold, Grid *energyGrid)
            : threshold(threshold), energyGrid(energyGrid) {}

    CrossSection::CrossSection(double threshold, Grid *energyGrid, std::ifstream &in)
            : threshold(threshold), energyGrid(energyGrid) {
        std::vector<double> rawEnergyVector, rawCrossSectionVector;

        Parse::rawCrossSectionFromStream(rawEnergyVector, rawCrossSectionVector, in);

        rawEnergyData = Vector::Map(rawEnergyVector.data(), rawEnergyVector.size());
        rawCrossSection = Vector::Map(rawCrossSectionVector.data(), rawCrossSectionVector.size());

        this->interpolate();
        this->energyGrid->updatedMaxEnergy1.addListener(&CrossSection::interpolate, this);
    }

    CrossSection::CrossSection(double threshold, Grid *energyGrid, Vector rawEnergyData, Vector rawCrossSection)
            : threshold(threshold), energyGrid(energyGrid), rawEnergyData(std::move(rawEnergyData)),
              rawCrossSection(std::move(rawCrossSection)) {

        this->interpolate();
        this->energyGrid->updatedMaxEnergy1.addListener(&CrossSection::interpolate, this);

    }

    void CrossSection::interpolate() {

        interpolate(energyGrid->getNodes(), *this);
    }

    void CrossSection::interpolate(const Vector &energies, Vector &result) {
        const uint32_t &gridSize = energies.size();

        result.resize(gridSize);

        uint32_t csIndex = 0;

        for (uint32_t gridIndex = 0; gridIndex < gridSize; ++gridIndex) {
            while (csIndex < rawCrossSection.size()
                   && rawEnergyData[csIndex] < energies[gridIndex])
                ++csIndex;

            if (csIndex >= rawCrossSection.size()) {
//                (*this)[gridIndex] = rawCrossSection.back().second;
                result[gridIndex] = 0.;
                continue;
            }

            if (csIndex == 0) {
                if (rawEnergyData[csIndex] == energies[gridIndex]) {
                    result[gridIndex] = rawCrossSection[csIndex];
                } else {
                    result[gridIndex] = 0.;
                }
                continue;
            }

            const double prevEnergy = rawEnergyData[csIndex - 1],
                    nextEnergy = rawEnergyData[csIndex],
                    prevCS = rawCrossSection[csIndex - 1],
                    nextCS = rawCrossSection[csIndex];

            const double alpha = (energies[gridIndex] - prevEnergy) / (nextEnergy - prevEnergy);

            result[gridIndex] = (1. - alpha) * prevCS + alpha * nextCS;
        }
    }

    Vector &CrossSection::raw() {
        return rawCrossSection;
    }

    Vector &CrossSection::energies() {
        return rawEnergyData;
    }
}