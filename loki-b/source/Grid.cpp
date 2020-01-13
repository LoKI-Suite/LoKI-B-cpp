//
// Created by daan on 2-5-19.
//

#include "../headers/Grid.h"
#include <Log.h>

namespace loki
{

Grid::Grid(const EnergyGridSetup &gridSetup)
    : cells(gridSetup.cellNumber), nodes(gridSetup.cellNumber - 1), 
      cellNumber(gridSetup.cellNumber), step(gridSetup.maxEnergy / gridSetup.cellNumber),
      minEedfDecay(gridSetup.smartGrid.minEedfDecay), maxEedfDecay(gridSetup.smartGrid.maxEedfDecay), 
      updateFactor(gridSetup.smartGrid.updateFactor), isSmart(updateFactor != 0 && minEedfDecay < maxEedfDecay)
{
    this->nodes = Vector::LinSpaced(cellNumber + 1, 0, gridSetup.maxEnergy);
    this->cells = Vector::LinSpaced(cellNumber, .5, cellNumber - .5) * step;

    Log<Message>::Notify(*this);
}

const Vector &Grid::getNodes() const
{
    return nodes;
}

const Vector &Grid::getCells() const
{
    return cells;
}

double Grid::getNode(uint32_t index) const
{
    return nodes[index];
}

double Grid::lastNode() const
{
    return nodes[cellNumber];
}

double Grid::getCell(uint32_t index) const
{
    return cells[index];
}

double Grid::lastCell() const
{
    return cells[cellNumber - 1];
}

void Grid::updateMaxEnergy(double value)
{
    step = value / cellNumber;
    nodes = Vector::LinSpaced(cellNumber + 1, 0, value);
    cells = Vector::LinSpaced(cellNumber, .5, cellNumber - .5) * step;

    updatedMaxEnergy1.emit();
    updatedMaxEnergy2.emit();
}

std::ostream &operator<<(std::ostream &os, const Grid &grid)
{
    return os << "Grid with " << grid.nodes.size() << " nodes, " << grid.cells.size()
              << " cells and step " << grid.step;
}
} // namespace loki