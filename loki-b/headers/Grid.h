//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_GRID_H
#define LOKI_CPP_GRID_H

#include "LinearAlgebra.h"
#include "Event.h"
#include "Setup.h"

namespace loki {
    class Grid {
        Vector cells, nodes;
        double step;

        // Smart grid
        bool isSmart = false;
        uint16_t minEedfDecay, maxEedfDecay;
        double updateFactor;

    public:
        uint32_t cellNumber;

        explicit Grid(const EnergyGridSetup &gridSetup);
        ~Grid() = default;

        Grid(const Grid &other) = delete;

        const Vector &getNodes() const;

        //Events
        event updatedMaxEnergy1,
              updatedMaxEnergy2;

        friend std::ostream& operator<<(std::ostream&, const Grid&);
    };
}


#endif //LOKI_CPP_GRID_H
