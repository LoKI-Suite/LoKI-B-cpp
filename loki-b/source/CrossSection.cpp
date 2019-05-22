//
// Created by daan on 15-5-19.
//

#include "CrossSection.h"
#include "Parse.h"

namespace loki {

    CrossSection::CrossSection(const double threshold, Grid *energyGrid) :
            threshold(threshold), energyGrid(energyGrid) {}

    CrossSection::CrossSection(double threshold, const Grid *energyGrid, std::ifstream &in)
            : threshold(threshold), energyGrid(energyGrid) {
        Parse::rawCrossSectionFromStream(rawCrossSection, in);
    }
}