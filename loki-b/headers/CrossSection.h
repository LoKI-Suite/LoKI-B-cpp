//
// Created by daan on 15-5-19.
//

#ifndef LOKI_CPP_CROSSSECTION_H
#define LOKI_CPP_CROSSSECTION_H

#include "LinearAlgebra.h"
#include "Grid.h"
#include "json.h"

#include <vector>
#include <fstream>

namespace loki {
    class CrossSection : public Vector {
//        std::vector<std::pair<double, double>> rawCrossSection;
        Vector rawEnergyData, rawCrossSection;
        Grid *energyGrid;

        const bool isElasticOrEffective;

    public:
        const double threshold;

        CrossSection(double threshold, Grid *energyGrid, bool isElasticOrEffective, const json_type& cnf);

        CrossSection(double threshold, Grid *energyGrid, bool isElasticOrEffective, std::ifstream &in);

        CrossSection(double threshold, Grid *energyGrid, bool isElasticOrEffective, Vector rawEnergyData,
                     Vector rawCrossSection);

        void interpolate();

        void interpolate(const Vector &energies, Vector &result);

        Vector &raw();

        Vector &energies();

        const Grid *getGrid();

//        CrossSection(double threshold, Grid *energyGrid, bool isElasticOrEffective);

    };
}

#endif //LOKI_CPP_CROSSSECTION_H
