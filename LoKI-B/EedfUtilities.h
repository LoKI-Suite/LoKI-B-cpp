/** \file
 *
 *  Interface of LoKI-B's linear algebra facilities.
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
 *  \date   29 July 2023
 */

#ifndef LOKI_CPP_EEDFUTILITIES_H
#define LOKI_CPP_EEDFUTILITIES_H

#include "LoKI-B/Exports.h"
#include "LoKI-B/Grid.h"

namespace loki {

/** Scale the elements of the edf such that the sum of sqrt(u_i)edf[i]du
 *  equals unity.
 */
lokib_export void normalizeEDF(Vector& edf, const Grid& grid);

/** Sample the prescribed EDF with shape factor \a g and temperature T_eV (in eV)
 *  on the cells of the provided \a grid and write the results in \a edf. Before
 *  returning, the resulting EDF is normalized by calling normalizeEDF(edf,grid).
 *
 *  The vector \a edf is resized to the number of grid cells by this function.
 *
 *  \note For g=1, a Maxwell EDF is obtained, for, g=2 a Druyvesteyn EDF.
 */
lokib_export void makePrescribedEDF(Vector& edf, const Grid& grid, double g, double T_eV);

/** See the overload of this function that accepts an edf parameter. The present
 *  function creates an edf vector, calls that overload to set the values and
 *  returns the result.
 */
lokib_export Vector makePrescribedEDF(const Grid& grid, double g, double T_eV);

} // namespace loki


#endif // LOKI_CPP_EEDFUTILITIES_H
