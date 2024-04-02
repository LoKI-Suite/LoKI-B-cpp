/** \file
 *
 *  Functions for querying environment variables.
 *  containers (per gas and for a mixture) of such collisions.
 *  This code is based on the PLASIMO code plutil/environment.h
 *  and accompanying .cpp file, see see https://plasimo.phys.tue.nl.
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
 *  \date   May 2013
 */

#ifndef LOKI_CPP_ENVIRONMENT_H
#define LOKI_CPP_ENVIRONMENT_H

#include <string>

namespace loki {

/** If the environment variable \a envvar exists, return its
 *  value as a std::string. If the variable is not set, and
 *  \a required is false, return an empty string is. If the
 *  variable is not set and \a required is true, throw a
 *  std::runtime_error.
 *
 *  \author Jan van Dijk
 *  \date   May 2013
 */
std::string getEnvironmentVariable(const std::string& envvar, bool required);

} // namespace loki

#endif // LOKI_CPP_ENVIRONMENT_H

