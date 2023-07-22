/** \file
 *
 *  Interface of the LoKI-B energy grid class.
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

/** The (energy) Grid spans an energy interval [0,u_max]. This region is
 *  divided into Nc cells of size du, where du = u_max/Nc. The position
 *  (energy) u_i of cell i is the energy at its centre and is given by
 *  (i+0.5)*du, where i is in the interval [0,Nc-1]. The cells are bounded by
 *  Nc+1 nodes at locations u_i^node = i*du, where i is in the interval [0,Nc].
 *
 *  The grid can be viewed schematically as
 *  \verbatim
        | x | x | ... | x | \endverbatim
 *
 *  where the | indicate nodes and the x cells.
 *
 *  The Grid class stores the energies in the nodes and cells in vectors.
 *  Constant references to these vectors can be retrieved with the members
 *  getNodes() and getCells(), respectively. In addition, members getNode
 *  and getCell allow the retrievel of the energy of a particlar node or
 *  cell for a given index. Members nCells() and uMax() return the number
 *  of cells and the energy of the last node.
 *
 *  Whereas the number of cells is fixed at construction time, the upper
 *  energy boundary can be modified later on by passing a new value in
 *  a call to the member function updateMaxEnergy. This functionality
 *  is used to implement the 'smart grid' functionality, in which the
 *  grid is adjusted to achieve a minimal or maximal dynamic range of the
 *  EEDF. The smart grid functionality is configured at construction time
 *  when a smartGrid section is present in the grid configuration. Access
 *  to the smart grid configuration is provided by the member function
 *  smartGrid(), which returns a pointer to the SmartGridParameters, or
 *  nullptr if no smart grid has been configured.
 *
 *  When the grid is changed during the simulation, various actions may
 *  need to be undertaken. This can be achieved by registering listeners to
 *  the event object updatedMaxEnergy. Those functions will be invoked at
 *  the end of a call to updateMaxEnergy().
 *
 *  \todo investigate the usage of the word node for a
 *        cell boundary, that sounds odd.
 *
 *  \author Daan Boer and Jan van Dijk (C++ version)
 *  \date   2. May 2019
 */
class Grid
{
public:

    /// The type of the vectors that hold the energy values
    using Vector = loki::Vector;
    /// The type of the index of elements of the energy vectors
    using Index = Vector::Index;

    // constructorion and destruction:

    Grid(unsigned nCells, double maxEnergy);
    explicit Grid(const EnergyGridSetup &gridSetup);
    /** Construct a Grid from the parameters "maxEnergy" (double)
     *  and "cellNumber" (unsigned) in the json object \a cnf.
     *  When an element "smartGrid" is present, a SmartGridParameters
     *  object will be created, see smartGrid().
     */
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

    /** Reconfigure the grid, using \a uMax as the new maximum energy.
     *  This function recalculates the energy values in the nodes and cells,
     *  then fires the event updatedMaxEnergy.
     */
    void updateMaxEnergy(double uMax);

    /** Event. This is mutable, so a new listener can be added to a constant
     *  grid object.
     */
    mutable Event<> updatedMaxEnergy;

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
    private:
        void checkConfiguration() const;
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
