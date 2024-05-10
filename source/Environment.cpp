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

#include "LoKI-B/Environment.h"
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace loki {

std::string getEnvironmentVariable(const std::string& envvar)
{
	const char* value = std::getenv(envvar.c_str());
	if (!value)
	{
		throw std::runtime_error( "Environment variable '"
			+ envvar + "' required, but not set.");
	}
	return value;
}

} // namespace loki
