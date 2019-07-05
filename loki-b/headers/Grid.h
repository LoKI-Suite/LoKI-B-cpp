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

    public:
        // Smart grid
        const uint16_t minEedfDecay, maxEedfDecay;
        const double updateFactor;
        const bool isSmart;

        uint32_t cellNumber;
        double step;

        explicit Grid(const EnergyGridSetup &gridSetup);
        ~Grid() = default;

        Grid(const Grid &other) = delete;

        const Vector &getNodes() const;

        const Vector &getCells() const;

        const double getNode(uint32_t index) const;

        const double lastNode() const;

        const double getCell(uint32_t index) const;

        const double lastCell() const;

        void updateMaxEnergy(double value);

        //Events
        event updatedMaxEnergy1,
              updatedMaxEnergy2;

        friend std::ostream& operator<<(std::ostream&, const Grid&);
    };
}


#endif //LOKI_CPP_GRID_H
