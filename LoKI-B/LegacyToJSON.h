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

/** Read properties from file \a fname. That file must be a LoKI-B
 *  database file. Percentage signs start comments, and lines that are
 *  non-empty after stripping the comments must contain key-value pairs,
 *  the values must be valid double values.
 *  The properties are stored as elements of a JSON object that contains
 *  the values for all keys and is returned by this function. A runtime
 *  error is thrown when the file cannot be read or contains bad data.
 */
json_type readLegacyGasPropertyFile(const std::filesystem::path& fname);

/** Parses legacy state property file \a fileName, stores the result in JSON
 *  format and returns the result. For the format of the JSON object, see the
 *  documentation of \c GasMixture::loadStateProperty.
 *
 *  After stripping LoKI-style comments, starting with a percentage sign
 *  the file must contain lines with two items, separated by whitespace,
 *  representing a state string and the initializer. The latter must be
 *  a numerical value or a function call. The latter takes the form of a
 *  string, followed by an argument list of the form "@arg1[,arg2]...",
 *  where the arguments are numerical values of parameter names. In order
 *  to support those, we can re-use the code from function patchStateProperty.
 */
lokib_export json_type readLegacyStatePropertyFile(const std::filesystem::path &fileName);

} // namespace loki

#endif // LOKI_CPP_LEGACYTOJSON_H
