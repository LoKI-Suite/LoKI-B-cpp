/** \file
 *
 *  Interface of LoKI-B EEDF-related utility functions.
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
 *  \author Daan Boer and Jan van Dijk (C++ version)
 *  \date   29 July 2023
 */

#ifndef LOKI_CPP_EEDFUTILITIES_H
#define LOKI_CPP_EEDFUTILITIES_H

#include "LoKI-B/Exports.h"
#include "LoKI-B/Grid.h"

namespace loki {

/** Calculate and return the mean energy in eV of particles that are distributed
 *  according an \a eedf that is defined on the cells of a \a grid. This is
 *  equal to \f$ \int_0^\infty u^{3/2}f(u)du \f$ and approximated by the sum
 *  of \f$ u_i^(3/2)f(u_i)(\Delta u)_i \f$ over the cells i.
 */
lokib_export double getMeanEnergy(const Vector& edf, const Grid& grid);

/** Scale the elements of the edf such that the sum of sqrt(u_i)edf[i](Delta u)_i
 *  equals unity.
 */
lokib_export void normalizeEDF(Vector& edf, const Grid& grid);

/** Sample the prescribed EDF with shape factor \a s and temperature T_eV (in
 *  eV) on the cells of the provided \a grid and write the results in \a edf.
 *  See \cite Nam2008 (eq. 3) for the version implemented here. This function
 *  also appears in the earlier work of \cite Amemiya1997 (eq. 13), in a
 *  slightly different form.
 *
 *  The vector \a edf is resized to the number of grid cells by this function.
 *  If \a normalize is true (the default), the resulting EDF is normalized by
 *  calling normalizeEDF(edf,grid). Note that the sampled values of the EEDF
 *  will deviate from the analytical EEDF function by a constant factor when
 *  normalization is enabled, whereas the computed approximation of the integral
 *  of \f$ \sqrt{u}f(u) \f$ over the grid's energy range will be unity.
 *
 *  \note For s=1 a Maxwell EDF is obtained, for s=2 a Druyvesteyn EDF.
 *  \note The mean energy of the resulting edf will not be exactly equal to
 *        (3/2)*T_eV because of discretisation errors. The normalization fixes
 *        the integral chance of finding a particle (==1), not the mean energy
 *        of the distribution.
 *
 *  \todo The risk of underflows should be investigated. (And diagnosed?)
 */
lokib_export void makePrescribedEDF(Vector& edf, const Grid& grid, double g, double T_eV, bool normalize=true);

/** See the overload of this function that accepts an edf parameter. The present
 *  function creates an edf vector, calls that overload to set the values and
 *  returns the result.
 */
lokib_export Vector makePrescribedEDF(const Grid& grid, double g, double T_eV, bool normalize=true);

} // namespace loki


#endif // LOKI_CPP_EEDFUTILITIES_H
