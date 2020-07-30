//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_GRID_H
#define LOKI_CPP_GRID_H

#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/Event.h"
#include "LoKI-B/Setup.h"
#include "LoKI-B/json.h"

namespace loki
{
class Grid
{
    Vector cells, nodes;

public:
    uint32_t cellNumber;
    double step;

    // Smart grid
    const uint16_t minEedfDecay, maxEedfDecay;
    const double updateFactor;
    const bool isSmart;

    explicit Grid(const EnergyGridSetup &gridSetup);
    explicit Grid(const json_type &cnf);
    ~Grid() = default;

    Grid(const Grid &other) = delete;

    const Vector &getNodes() const;

    const Vector &getCells() const;

    double getNode(uint32_t index) const;

    double lastNode() const;

    double getCell(uint32_t index) const;

    double lastCell() const;

    void updateMaxEnergy(double value);

    //Events
    event updatedMaxEnergy1,
        updatedMaxEnergy2;

    friend std::ostream &operator<<(std::ostream &, const Grid &);
};
} // namespace loki

#endif //LOKI_CPP_GRID_H
