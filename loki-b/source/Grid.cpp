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

    std::ostream &operator<<(std::ostream &os, const Grid &grid) {
        return os << "Grid has " << grid.nodes.size() << " nodes, " << grid.cells.size()
                  << " cells and step " << grid.step;
    }
}