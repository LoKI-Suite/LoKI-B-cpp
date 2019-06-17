//
// Created by daan on 2-5-19.
//

#include "../headers/Grid.h"
#include <Log.h>

namespace loki {


    Grid::Grid(const EnergyGridSetup &gridSetup)
        : cellNumber(gridSetup.cellNumber), cells(gridSetup.cellNumber),
            nodes(gridSetup.cellNumber-1) {

        this->step = gridSetup.maxEnergy / cellNumber;

        // TODO: Vector::LinSpaced is EIGEN specific.
        this->nodes = Vector::LinSpaced(cellNumber+1, 0, gridSetup.maxEnergy);
        this->cells = Vector::LinSpaced(cellNumber, .5, cellNumber-.5) * step;

        // Smart grid
        this->updateFactor = gridSetup.smartGrid.updateFactor;
        this->minEedfDecay = gridSetup.smartGrid.minEedfDecay;
        this->maxEedfDecay = gridSetup.smartGrid.maxEedfDecay;

        if (updateFactor != 0 && minEedfDecay < maxEedfDecay) {

            this->isSmart = true;
        }

        Log<Message>::Notify(*this);
    }

    const Vector &Grid::getNodes() const {
        return nodes;
    }

    const Vector &Grid::getCells() const {
        return cells;
    }

    const double Grid::getNode(uint32_t index) const {
        return nodes[index];
    }

    const double Grid::getCell(uint32_t index) const {
        return cells[index];
    }

    std::ostream &operator<<(std::ostream &os, const Grid &grid) {
        return os << "Grid with " << grid.nodes.size() << " nodes, " << grid.cells.size()
                  << " cells and step " << grid.step;
    }
}