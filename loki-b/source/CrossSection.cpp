//
// Created by daan on 15-5-19.
//

#include "CrossSection.h"
#include "Parse.h"
#include "Log.h"

namespace loki {

    CrossSection::CrossSection(const double threshold, Grid *energyGrid) :
            threshold(threshold), energyGrid(energyGrid) {}

    CrossSection::CrossSection(double threshold, Grid *energyGrid, std::ifstream &in)
            : threshold(threshold), energyGrid(energyGrid) {
        Parse::rawCrossSectionFromStream(rawCrossSection, in);

        this->interpolate();
        this->energyGrid->updatedMaxEnergy1.addListener(&CrossSection::interpolate, this);
    }

    void CrossSection::interpolate() {
        const auto &nodes = energyGrid->getNodes();
        const uint32_t &gridSize = nodes.size();

        this->resize(gridSize);

        uint32_t csIndex = 0;

        for (uint32_t gridIndex = 0; gridIndex < gridSize; ++gridIndex) {
            while(rawCrossSection[csIndex].first < nodes[gridIndex]
                  && csIndex < rawCrossSection.size()) ++csIndex;

            if (csIndex >= rawCrossSection.size()) {
//                (*this)[gridIndex] = rawCrossSection.back().second;
                (*this)[gridIndex] = 0.;
                continue;
            }

            if (csIndex == 0) {
                if (rawCrossSection[csIndex].first == nodes[gridIndex]) {
                    (*this)[gridIndex] = rawCrossSection[csIndex].second;
                }
                else {
                    (*this)[gridIndex] = 0.;
                }
                continue;
            }

            const auto &prev = rawCrossSection[csIndex-1],
                    &next = rawCrossSection[csIndex];

            const double alpha = (nodes[gridIndex] - prev.first) / (next.first - prev.first);

            (*this)[gridIndex] = (1. - alpha) * prev.second + alpha * next.second;
        }
    }
}