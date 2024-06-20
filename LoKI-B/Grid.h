/** \file
 *
 *  Interface of the LoKI-B energy grid class.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2024 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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
 *  \author Daan Boer, Jan van Dijk and Jop Hendrikx (C++ version)
 *  \date   2 May 2019 (initial version)
 */

#ifndef LOKI_CPP_GRID_H
#define LOKI_CPP_GRID_H

#include "LoKI-B/Event.h"
#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/json.h"
#include <memory>
#include <stdexcept>

namespace loki
{

/** \brief The (energy) Grid spans an energy interval [0,u_max].
 *
 *  This region is divided into Nc cells that are bounded by Nf=Nc+1 faces at
 *  energies uFaces_f. The grid can be viewed schematically as
 *
 *  \verbatim
        0   1     2     Nf-2    Nf-1
        | x |  x  | ...   |   x   |
          0    1            Nc-1 \endverbatim
 *
 *  Here the | indicate faces and the x cells. The indices above the face
 *  locations are the indices f of the faces, 0<=f<Nf, the indices below
 *  below the cell locations represent those of the cells, 0<=c<Nc.
 *
 *  From the picture it is clear that cell c is bounded by the faces at indices
 *  c and c+1. The energy uCells[c] of this cell is taken to be that in in the
 *  middle of the cell, uCell[c] := (uFaces[c] + uFaces[c+1])/2.
 *
 *  In general, the grid is not uniform. This means that the width of cell c,
 *  which is given by duCell[c] = uFaces[c+1]-uFaces[c], depends on c.
 *
 *  The 'width' duFaces[f] of an internal interface, 1<=f<Nf-1, is defined as
 *  the energy distance between the adjacent cells. For boundary faces, it is
 *  defined as the distance from the boundary to the adjacent cell). We then
 *  have: \verbatim
    f=0:       duFaces[f] = uCells[f]
    1<=f<Nf-1: duFaces[f] = uCells[f] - uCells[f-1]
    f=Nf-1:    duFaces[f] = u_max - uCells[f] \endverbatim
 *
 *  In the special case of a uniform energy grid, duCells[c] = u_Max/Nc := dU,
 *  and duFaces[f] = dU for internal faces, dU/2 for boundary faces.
 *
 *  A grid object can be constructed in various ways. A constructor overload
 *  that accepts arguments nCells and u_max will set up a uniform grid with the
 *  given number of cells and energy range [0,u_max]. Another overload
 *  creates a grid from a collection of sorted normalized node positions,
 *  which must span the range [0,1], and u_max. This will in general result
 *  in a non-uniform mesh. Another constructor overloads create a grid object
 *  from the parameters that are specified in a json_type object.
 *
 *  The Grid class stores the energies of the faces and cells in vectors.
 *  Constant references to these vectors can be retrieved with the members
 *  getNodes() and getCells(), respectively. In addition, overloads of members
 *  getNode and getCell that return energy for a face or cell at a given index
 *  are available. Similarly, members duNodes and duCells give access to
 *  the widths of the faces and cells. Members nCells() and uMax() return the
 *  number of cells and the energy of the last node. Member isUniform()
 *  returns a boolean that can be inspected to see whether the grid is uniform:
 *  this may allow some grid-optimizations to be optimized.
 *
 *  Whereas the number of cells and the relative distribution of the face
 *  energies are fixed at construction time, the upper energy boundary can be
 *  modified later on by passing the new value as an argument to member function
 *  updateMaxEnergy.
 *
 *  Finally, a grid object manages the parameters that are used by the 'Smart
 *  Grid' functionality. These can be obtained with the member function
 *  smartGrid(), which returns a pointer to an object of the nested class type
 *  SmartGridParameters, or a nullptr when no smart grid support was configured.
 *  The smart grid functionality is configured at construction time when a
 *  smartGrid section is present in the grid configuration.
 *  See the SmartGridParameters documentation for details.
 *
 *  When the grid is changed during the simulation, various actions may
 *  need to be undertaken. This can be achieved by registering listeners to
 *  the event object updatedMaxEnergy. Those functions will be invoked at
 *  the end of a call to updateMaxEnergy().
 *
 *  \todo investigate the usage of the word node for a
 *        cell boundary, that sounds odd.
 *
 *  \author Daan Boer, Jan van Dijk and Jop Hendrikx (C++ version)
 *  \date   2 May 2019 (initial version)
 */

class EedfMixture;

class Grid
{
  public:
    /// The type of the vectors that hold the energy values
    using Vector = loki::Vector;
    /// The type of the index of elements of the energy vectors
    using Index = Vector::Index;

    // Forward declaration
    class SmartGridParameters;

    // construction and destruction:

    Grid(unsigned nCells, double maxEnergy, SmartGridParameters *smartGridParameters = nullptr);
    /** Construct a nonuniform Grid from the parameters "nodeDistribution"
     * (Vector) and "maxEnergy" (double).
     * "nodeDistribution" contains doubles in the range [0-1] which indicates
     * the normalized distribution of cell faces (nodes) over the energy
     * domain. "maxEnergy" is the maximum energy of the simulation domain.
     * The physical energy values of the nodes can be calculated by multiplying
     * the maxEnergy and the nodeDistribution function.
     */
    Grid(const Vector &nodeDistribution, double maxEnergy);
    /** Construct a Grid from the parameters "maxEnergy" (double)
     *  and "cellNumber" (unsigned) in the json object \a cnf.
     *  When an element "smartGrid" is present, a SmartGridParameters
     *  object will be created, see smartGrid().
     */
    explicit Grid(const json_type &cnf);
    Grid(const Vector &nodeDistribution, double maxEnergy, bool isUniform);
    /// Grids canot be copied.
    Grid(const Grid &other) = delete;

    static Grid fromConfig(const json_type &cnf);
    /// The default destructor is used for Grid objects
    ~Grid() = default;

    // accessors:

    double uMax() const
    {
        return m_nodes[nCells()];
    }
    Index nCells() const
    {
        return m_nCells;
    }
    double du() const
    {
        if (!m_isUniform)
        {
            throw std::runtime_error("Grid::du() is not available for non-uniform meshes.");
        }
        return m_du;
    }
    const Vector &getNodes() const
    {
        return m_nodes;
    }
    const Vector &getCells() const
    {
        return m_cells;
    }
    /** Vector of 'face widths': the distances between a face' adjacent cells.
     *  For boundary faces, the distance between the energy at the boundary
     *  and that of the adjacent internal cell is used.
     */
    const Vector &duNodes() const
    {
        return m_duNodes;
    }
    /// Vector of cell widths (one for each cell).
    const Vector &duCells() const
    {
        return m_duCells;
    }
    double getNode(Index index) const
    {
        return m_nodes[index];
    }
    double getCell(Index index) const
    {
        return m_cells[index];
    }

    bool isUniform() const
    {
        return m_isUniform;
    }

    double duNode(Index index) const
    {
        return m_duNodes[index];
    }

    double duCell(Index index) const
    {
        return m_duCells[index];
    }

    /** Reconfigure the grid, using \a uMax as the new maximum energy.
     *  This function recalculates the energy values in the nodes and cells,
     *  then fires the event updatedMaxEnergy.
     */
    void updateMaxEnergy(double uMax);

    void updateMaxEnergyNonuniform(double uMax, const EedfMixture &mixture);

    /** Event. This is mutable, so a new listener can be added to a constant
     *  grid object.
     */
    mutable Event<> updatedMaxEnergy;

    /** A SmartGridParameters object is used to implement the 'smart grid'
     *  functionality. This means that the u_max of the  grid is adjusted to
     *  achieve aan EEDF with a dynamic range within a user-specified interval.
     *  The smart grid functionality is configured at construction time
     *  when a smartGrid section is present in the grid configuration.
     */
    class SmartGridParameters
    {
      public:
        SmartGridParameters(const json_type &cnf);
        /// \todo We don't we use a more straighforward type here, like unsigned?
        const uint16_t minEedfDecay;
        const uint16_t maxEedfDecay;
        const double updateFactor;

      private:
        void checkConfiguration() const;
    };

    const SmartGridParameters *smartGrid() const
    {
        return m_smartGrid.get();
    }

  private:
    const Index m_nCells;
    /// the energy spacing (only used for a uniform grid)
    double m_du;
    /// energies at the faces (nodes)
    Vector m_nodes;
    /// energies at the cell centers
    Vector m_cells;
    std::unique_ptr<const SmartGridParameters> m_smartGrid;
    /// indicates if the grid is uniform
    bool m_isUniform;
    /** Vector of 'face widths': the distances between a face' adjacent cells.
     *  For boundary faces, the distance between the energy at the boundary
     *  and that of the adjacent internal cell is used.
     */
    Vector m_duNodes;
    /// Vector of cell widths (one for each cell).
    Vector m_duCells;
};

} // namespace loki

#endif // LOKI_CPP_GRID_H
