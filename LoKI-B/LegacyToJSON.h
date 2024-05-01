/** \file
 *
 *  Conversion of legacy LoKI-B input files to LoKI-B JSON format.
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
 *  \author Jan van Dijk
 *  \date   30 April 2024
 */

#ifndef LOKI_CPP_LEGACYTOJSON_H
#define LOKI_CPP_LEGACYTOJSON_H

#include "LoKI-B/json.h"
#include <iostream>
#include <filesystem>

namespace loki
{

/** This function makes a copy of a legacy LoKI-B input file, converted to a
 *  JSON object using offSideToJSON, applies the necessary patches to make the
 *  result a valid LoKI-B JSON input file and returns the result. See section
 *  5 of \cite Manual_2_2_0 for details about the legacy file format.
 *
 *  \todo Document the required, JSON format, provide a reference here.
 *
 *  \sa offSideToJSON
 *  \author Jan van Dijk
 *  \date   30 April 2024
 */
json_type legacyToJSON(const json_type& legacy);

/** Reads an 'off-side rule' file from stream \a is, converts it to a JSON
 *  object using function offSideToJSON, process it using the json_type overload
 *  of this function and return the result. 
 *
 *  \sa offSideToJSON
 *  \author Jan van Dijk
 *  \date   30 April 2024
 */
json_type legacyToJSON(std::istream& is);

/** Creates and input file stream for file \a fname and returns the result of
 *  calling the stream-overload of this function. A runtime_error is thrown if
 *  file \a fname could not be opened.
 *
 *  \sa offSideToJSON
 *  \author Jan van Dijk
 *  \date   30 April 2024
 */
json_type legacyToJSON(const std::filesystem::path& fname);

} // namespace loki



#endif // LOKI_CPP_LEGACYTOJSON_H

