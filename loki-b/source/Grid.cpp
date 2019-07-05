//
// Created by daan on 2-5-19.
//

#include "../headers/Grid.h"
#include <Log.h>

namespace loki {


    Grid::Grid(const EnergyGridSetup &gridSetup)
            : cellNumber(gridSetup.cellNumber), cells(gridSetup.cellNumber),
              nodes(gridSetup.cellNumber - 1), updateFactor(gridSetup.smartGrid.updateFactor),
              minEedfDecay(gridSetup.smartGrid.minEedfDecay), maxEedfDecay(gridSetup.smartGrid.maxEedfDecay),
              isSmart(updateFactor != 0 && minEedfDecay < maxEedfDecay) {

        this->step = gridSetup.maxEnergy / cellNumber;

        // TODO: Vector::LinSpaced is EIGEN specific.
        this->nodes = Vector::LinSpaced(cellNumber + 1, 0, gridSetup.maxEnergy);
        this->cells = Vector::LinSpaced(cellNumber, .5, cellNumber - .5) * step;

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

    const double Grid::lastNode() const {
        return nodes[cellNumber];
    }

    const double Grid::getCell(uint32_t index) const {
        return cells[index];
    }

    const double Grid::lastCell() const {
        return cells[cellNumber - 1];
    }

    void Grid::updateMaxEnergy(double value) {
        step = value / cellNumber;
        nodes = Vector::LinSpaced(cellNumber + 1, 0, value);
        cells = Vector::LinSpaced(cellNumber, .5, cellNumber - .5) * step;

        updatedMaxEnergy1.emit();
        updatedMaxEnergy2.emit();
    }

    std::ostream &operator<<(std::ostream &os, const Grid &grid) {
        return os << "Grid with " << grid.nodes.size() << " nodes, " << grid.cells.size()
                  << " cells and step " << grid.step;
    }
}