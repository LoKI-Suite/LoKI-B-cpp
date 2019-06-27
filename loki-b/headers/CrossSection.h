//
// Created by daan on 15-5-19.
//

#ifndef LOKI_CPP_CROSSSECTION_H
#define LOKI_CPP_CROSSSECTION_H

#include "LinearAlgebra.h"
#include "Grid.h"

#include <vector>
#include <fstream>

namespace loki {
    class CrossSection : public Vector {
//        std::vector<std::pair<double, double>> rawCrossSection;
        Vector rawEnergyData, rawCrossSection;
        Grid *energyGrid;

    public:
        const double threshold;

        CrossSection(double threshold, Grid *energyGrid, std::ifstream &in);
        CrossSection(double threshold, Grid *energyGrid, Vector rawEnergyData, Vector rawCrossSection);

        void interpolate();
        void interpolate(const Vector &energies, Vector &result);

        Vector &raw();

        Vector &energies();

        const Grid *getGrid();

        CrossSection(double threshold, Grid *energyGrid);

    };
}

#endif //LOKI_CPP_CROSSSECTION_H
