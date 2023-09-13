/** \file
 *
 *  Implementation of the LoKI-B energy grid class.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2020 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
 *  M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
 *  L. L. Alves
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  \author Daan Boer and Jan van Dijk (C++ version)
 *  \date   2. May 2019
 */

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

Grid::Grid(unsigned nCells, double maxEnergy)
    : m_nCells(nCells),
      m_du(maxEnergy/m_nCells),
      m_nodes(Vector::LinSpaced(m_nCells + 1, 0, maxEnergy)),
      m_cells(Vector::LinSpaced(m_nCells, .5, m_nCells - .5) * m_du),  
      m_isUniform(true)
{
}

Grid::Grid(Vector nodeDistribution, double maxEnergy)
   : m_nCells(nodeDistribution.size()),
   m_du(0),
   m_nodes(nodeDistribution*maxEnergy),
   m_cells(.5*(m_nodes.tail(m_nCells - 1) + m_nodes.head(m_nCells - 1))),
   m_isUniform(false),
   m_duNodes(m_nodes.tail(m_nCells - 1) - m_nodes.head(m_nCells - 1)),
   m_duCells(m_cells.tail(m_nCells - 1) - m_cells.head(m_nCells - 1))
{
}

Grid::Grid(const EnergyGridSetup &gridSetup)
    : Grid(gridSetup.cellNumber,gridSetup.maxEnergy)
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
    : Grid(cnf.at("cellNumber").get<unsigned>(),cnf.at("maxEnergy").get<double>())
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
