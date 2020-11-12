//
// Created by daan on 2-5-19.
//

#include "LoKI-B/Grid.h"
#include "LoKI-B/Log.h"

namespace loki
{

Grid::SmartGridParameters::SmartGridParameters(const SmartGridSetup& smartGrid)
    : minEedfDecay(smartGrid.minEedfDecay),
      maxEedfDecay(smartGrid.maxEedfDecay), updateFactor(smartGrid.updateFactor),
      isSmart(updateFactor != 0 && minEedfDecay < maxEedfDecay)
{
    /// \todo Check parameters and throw in case of an error
}

Grid::SmartGridParameters::SmartGridParameters(const json_type& cnf)
    : minEedfDecay(cnf.at("minEedfDecay").get<unsigned>()),
      maxEedfDecay(cnf.at("maxEedfDecay").get<unsigned>()),
      updateFactor(cnf.at("updateFactor").get<double>()),
      isSmart(updateFactor != 0 && minEedfDecay < maxEedfDecay)
{
    /// \todo Check parameters and throw in case of an error
}

Grid::Grid(const EnergyGridSetup &gridSetup)
    : m_nCells(gridSetup.cellNumber),
      m_du(gridSetup.maxEnergy/m_nCells),
      m_nodes(Vector::LinSpaced(m_nCells + 1, 0, gridSetup.maxEnergy)),
      m_cells(Vector::LinSpaced(m_nCells, .5, m_nCells - .5) * m_du)
{
    /// \todo Set up m_smartGrid only if a smartGrid section has been specified
    m_smartGrid.reset(new SmartGridParameters(gridSetup.smartGrid));
    /** \todo Avoid creation of the structure in the first place
     *        Require valid parameters when any are specified.
     */
    if (!m_smartGrid->isSmart)
    {
        m_smartGrid.reset(nullptr);
    }
    Log<Message>::Notify("Created the energy grid.");
}

/// \todo Use get<T>() everywhere, do not rely on implicit conversion
/// \todo See if the elements for smartGrid support can be bundled in some way
Grid::Grid(const json_type &cnf)
    : m_nCells(cnf.at("cellNumber")),
      m_du(cnf.at("maxEnergy").get<double>()/m_nCells),
      m_nodes(Vector::LinSpaced(m_nCells + 1, 0, cnf.at("maxEnergy"))),
      m_cells(Vector::LinSpaced(m_nCells, .5, m_nCells - .5) * m_du)
{
    if (cnf.contains("smartGrid"))
    {
        m_smartGrid.reset(new SmartGridParameters(cnf.at("smartGrid")));
    }
    /** \todo Avoid creation of the structure in the first place
     *        Require valid parameters when any are specified.
     */
    if (m_smartGrid.get() && !m_smartGrid->isSmart)
    {
        m_smartGrid.reset(nullptr);
    }
    Log<Message>::Notify("Created the energy grid.");
}

void Grid::updateMaxEnergy(double uMax)
{
    m_du = uMax / nCells();
    m_nodes = Vector::LinSpaced(nCells() + 1, 0, uMax);
    m_cells = Vector::LinSpaced(nCells(), .5, nCells() - .5) * m_du;

    updatedMaxEnergy1.emit();
    updatedMaxEnergy2.emit();
}

} // namespace loki
