//
// Created by daan on 2-5-19.
//

#include "LoKI-B/Grid.h"
#include "LoKI-B/Log.h"

namespace loki
{

Grid::SmartGridParameters::SmartGridParameters(const SmartGridSetup& smartGrid)
    : minEedfDecay(smartGrid.minEedfDecay),
      maxEedfDecay(smartGrid.maxEedfDecay), updateFactor(smartGrid.updateFactor)
{
    checkConfiguration();
}

Grid::SmartGridParameters::SmartGridParameters(const json_type& cnf)
    : minEedfDecay(cnf.at("minEedfDecay").get<unsigned>()),
      maxEedfDecay(cnf.at("maxEedfDecay").get<unsigned>()),
      updateFactor(cnf.at("updateFactor").get<double>())
{
    checkConfiguration();
}

void Grid::SmartGridParameters::checkConfiguration() const
{
    if (updateFactor<=0)
    {
        Log<Message>::Error("Smart grid: update factor must be positive.");
    }
    if (minEedfDecay > maxEedfDecay)
    {
        Log<Message>::Error("Smart grid: minEedfDecay exceeds maxEedfDecay.");
    }
}

Grid::Grid(const EnergyGridSetup &gridSetup)
    : m_nCells(gridSetup.cellNumber),
      m_du(gridSetup.maxEnergy/m_nCells),
      m_nodes(Vector::LinSpaced(m_nCells + 1, 0, gridSetup.maxEnergy)),
      m_cells(Vector::LinSpaced(m_nCells, .5, m_nCells - .5) * m_du)
{
    const SmartGridSetup& sg = gridSetup.smartGrid;
    /* Try to set up the smartGrid only if any of the configuration
     * parameters does not have the default value (see Setup.h).
     * We use that to test that a configuration has been read from file.
     */
    if (sg.minEedfDecay!=0 || sg.maxEedfDecay!=0 || sg.updateFactor !=0.0)
    {
        m_smartGrid.reset(new SmartGridParameters(gridSetup.smartGrid));
    }
    Log<Message>::Notify("Created the energy grid.");
}

Grid::Grid(const json_type &cnf)
    : m_nCells(cnf.at("cellNumber").get<unsigned>()),
      m_du(cnf.at("maxEnergy").get<double>()/m_nCells),
      m_nodes(Vector::LinSpaced(m_nCells + 1, 0, cnf.at("maxEnergy"))),
      m_cells(Vector::LinSpaced(m_nCells, .5, m_nCells - .5) * m_du)
{
    if (cnf.contains("smartGrid"))
    {
        m_smartGrid.reset(new SmartGridParameters(cnf.at("smartGrid")));
    }
    Log<Message>::Notify("Created the energy grid.");
}

void Grid::updateMaxEnergy(double uMax)
{
    m_du = uMax / nCells();
    m_nodes = Vector::LinSpaced(nCells() + 1, 0, uMax);
    m_cells = Vector::LinSpaced(nCells(), .5, nCells() - .5) * m_du;

    updatedMaxEnergy.emit();
}

} // namespace loki
