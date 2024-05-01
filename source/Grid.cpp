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

Grid::Grid(unsigned nCells, double maxEnergy, SmartGridParameters *smartGridParameters)
   : Grid(Vector::LinSpaced(nCells + 1, 0.0, 1.0), maxEnergy, true)
{
    if (smartGridParameters != nullptr) {
        m_smartGrid.reset(smartGridParameters);
    }
}

Grid::Grid(const Vector& nodeDistribution, double maxEnergy)
   : Grid(nodeDistribution,maxEnergy,false)
{
}

Grid::Grid(const Vector& nodeDistribution, double maxEnergy, bool isUniform)
 : m_nCells(nodeDistribution.size()-1),
   m_du(isUniform ? maxEnergy/(nodeDistribution.size()-1) : 0.0),
   m_nodes(nodeDistribution*maxEnergy),
   m_cells(.5*(m_nodes.tail(m_nodes.size() - 1) + m_nodes.head(m_nodes.size() - 1))),
   m_isUniform(isUniform),
   m_duNodes(m_nodes.size()), // note: values are set in the constructor body
   m_duCells(m_nodes.tail(m_nodes.size() - 1) - m_nodes.head(m_nodes.size() - 1))
{
    /* the first node (==face) size is the energy difference between the
     * grid boundary and the adjacent internal cell.
     */
    m_duNodes[0] = m_cells[0] - 0.;
    m_duNodes[m_nodes.size()-1] = maxEnergy - m_cells[m_cells.size()-1];
    // for internal nodes, du is the difference between the energies of the adjacent cells.
    m_duNodes.segment(1,m_nodes.size()-2) = m_cells.tail(m_nCells - 1) - m_cells.head(m_nCells - 1);
#if 0
	std::cout << "nCells = " << m_nCells << std::endl;
	std::cout << "faces = " << m_nodes << std::endl;
	std::cout << "cells = " << m_cells << std::endl;
#endif
}

Grid Grid::fromConfig(const json_type &cnf)
{
    if (cnf.contains("nonuniformGrid"))
    {   
        auto nodes = cnf.at("nonuniformGrid").at("nodeDistribution").get<std::vector<double>>();
        Vector nodeDistribution = Eigen::Map<Vector, Eigen::Unaligned>(nodes.data(), nodes.size());
        return Grid(nodeDistribution, cnf.at("nonuniformGrid").at("maxEnergy").get<double>());
    } else
    {
        if (cnf.contains("smartGrid")) {
            return Grid(
                cnf.at("cellNumber").get<unsigned>(),
                cnf.at("maxEnergy").get<double>(), 
                new SmartGridParameters(cnf.at("smartGrid"))
            );
        }
        return Grid(cnf.at("cellNumber").get<unsigned>(),cnf.at("maxEnergy").get<double>());
    }
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
    const double uMaxRatio = uMax / this->uMax();
    m_du *= uMaxRatio;
    m_nodes *= uMaxRatio;
    m_cells *= uMaxRatio;
    m_duNodes *= uMaxRatio;
    m_duCells *= uMaxRatio;

    updatedMaxEnergy.emit();
}

} // namespace loki
