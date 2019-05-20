//
// Created by daan on 15-5-19.
//

#ifndef LOKI_CPP_CROSSSECTION_H
#define LOKI_CPP_CROSSSECTION_H

#include "LinearAlgebra.h"

#include <vector>

namespace loki {
    class CrossSection : public Vector {
        std::vector<std::pair<double, double>> rawCrossSection;

    public:
        const double threshold;

        explicit CrossSection(double threshold);

    };
}

#endif //LOKI_CPP_CROSSSECTION_H
