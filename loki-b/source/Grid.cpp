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
    : cellNumber(gridSetup.cellNumber),
      step(gridSetup.maxEnergy / gridSetup.cellNumber),
      m_nodes(Vector::LinSpaced(cellNumber + 1, 0, gridSetup.maxEnergy)),
      m_cells(Vector::LinSpaced(cellNumber, .5, cellNumber - .5) * step)
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
    : cellNumber(cnf.at("cellNumber")),
      step(cnf.at("maxEnergy").get<double>() / cnf.at("cellNumber").get<unsigned>()),
      m_nodes(Vector::LinSpaced(cellNumber + 1, 0, cnf.at("maxEnergy"))),
      m_cells(Vector::LinSpaced(cellNumber, .5, cellNumber - .5) * step)
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
    step = uMax / cellNumber;
    m_nodes = Vector::LinSpaced(cellNumber + 1, 0, uMax);
    m_cells = Vector::LinSpaced(cellNumber, .5, cellNumber - .5) * step;

    updatedMaxEnergy1.emit();
    updatedMaxEnergy2.emit();
}

} // namespace loki
