/** \file
 *
 *  The LoKI-B energy grid class.
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

#ifndef LOKI_CPP_GRID_H
#define LOKI_CPP_GRID_H

#include "LoKI-B/Event.h"
#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/Setup.h"
#include "LoKI-B/json.h"
#include <memory>

namespace loki
{

/** The (energy) Grid spans an energy interval [0,u_max].
 *  This region is divided into Nc cells of size du, where
 *  du = u_max/Nc. The positions (energies) u_i of cell i
 *  is the central energy of cell i and is given by
 *  (i+0.5)*du, where i is in the interval [0,Nc-1].
 *  The cells are bounded by Nc+1 nodes at locations
 *  u_i^node = i*du, where i is in the interval [0,Nc];
 *
 *  Below: | indicates a node, x a cell. The index values
 *  of the first and last nodes and cells are also indicated.
 *
 *  node index: 0                 Nc+1
 *              v                 v
 *              | x | x | ... | x |
 *                ^             ^
 *  cell index:   0             Nc
 *
 *  \todo investigate the usage of the word node for a
 *        cell boundary, that sounds odd.
 */
class Grid
{
  public:

    using Vector = loki::Vector;
    /** \todo IMHO, Index should be defined as Vector::Index for
     *        reasons of consistency; this also requires
     *        changes in other places where uint32_t is used
     *        at present.
     */
    using Index = uint32_t;
    //using Index = Vector::Index;

    // constructorion and desrtuction:

    explicit Grid(const EnergyGridSetup &gridSetup);
    explicit Grid(const json_type &cnf);
    /// The default constructor is used for Grid objects
    ~Grid() = default;
    /// Grids canot be copied.
    Grid(const Grid &other) = delete;

    // accessors:

    double uMax() const { return m_nodes[nCells()]; }
    Index nCells() const { return m_nCells; }
    double du() const { return m_du; }
    const Vector &getNodes() const { return m_nodes; }
    const Vector &getCells() const { return m_cells; }
    double getNode(Index index) const { return m_nodes[index]; }
    double getCell(Index index) const { return m_cells[index]; }

    /** Reconfigure the grid, using \a uMax is the new maximum energy.
     *  This function recalculates the energy values in the nodes and cells,
     *  then fires the events updatedMaxEnergy1 and updatedMaxEnergy2,
     *  in that order.
     */
    void updateMaxEnergy(double uMax);

    /** Events
     *  \todo Why do we need two? Is it possible to let the names of the
     *        objects reflect the reason (if we really need two).
     */
    Event<> updatedMaxEnergy1, updatedMaxEnergy2;

    /** Smart grid parameters.
     */
    class SmartGridParameters
    {
    public:
        SmartGridParameters(const SmartGridSetup& smartGrid);
        SmartGridParameters(const json_type& cnf);
        /// \todo We don't we use a more straighforward type here, like unsigned?
        const uint16_t minEedfDecay;
        const uint16_t maxEedfDecay;
        const double updateFactor;
        /// \todo Obsolete. when !isSmart, we will simply not create the smartGrid ptr.
        const bool isSmart;
    };
    const SmartGridParameters* smartGrid() const { return m_smartGrid.get(); }

private:
    const Index m_nCells;
    // the energy spacing
    double m_du;
    /// energies at the nodes (cell boundaries)
    Vector m_nodes;
    /// energies at the cell centers
    Vector m_cells;
    std::unique_ptr<const SmartGridParameters> m_smartGrid;
};

} // namespace loki

#endif // LOKI_CPP_GRID_H
