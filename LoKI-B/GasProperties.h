/** \file
 *
 *  Code for reading gas properties and storing the results in a JSON object.
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
 *  \date   2. May 2024
 */

#ifndef LOKI_CPP_GASPROPERTIES_H
#define LOKI_CPP_GASPROPERTIES_H

#include "LoKI-B/json.h"
#include <filesystem>

namespace loki {

/** This class manages the information that is contained in a LoKI-B
 *  'gasProperties' section. It stores the data in a private json object.
 *  The elements' keys are the property names, the values are JSON objects.
 *  The latter's elements' keys are the names of the gases for which the
 *  property is available, the values are the values of that property.
 *
 *  \author Jan van Dijk
 *  \date   2. May 2024
 */
class GasProperties
{
public:
    /** Set up a default GasProperties object. This adds one property: the
     *  "mass" of "e", the electron and is set to loki::Constant::electronMass
     *  (See LoKI-B/Constant.h).
     */
    GasProperties();
    /** Read the properties that are specified in LoKI-B gasProperties section
     *  \a pnode. The \a basePath is the file relative to which included files
     *  will be resolved (e.g. Databases/masses.txt). The \a pnode must be a
     *  JSON object. The key of each element is used as property name, the
     *  argument can be:
     *  1) a file name, as in "mass": "Databases/masses.txt"
     *  2) an array of key-double pairs. As an example, the fractions of gases
     *     N2 and O2 can be configured by: "fraction": { "N2": 0.8, "O2": 0.2 }.
     */
    GasProperties(const std::filesystem::path &basePath, const json_type& pnode);
    /** Return true if this object has a collection (possibly empty) of values
     *  of property \a propertyName.
     */
    bool has(const std::string& propertyName) const;
    /** Returns true if property \a propertyName is available for item \a key.
     */
    bool has(const std::string& propertyName, const std::string& key) const;
    /** Get property \a propertyName for item \a key. If the property cannot be
     *  found, a runtime_exception is thrown.
     */
    double get(const std::string& propertyName, const std::string& key) const;
    /** Get property \a propertyName for item \a key. If the property cannot be
     *  found, return alternative value \a alt instead. When \a warn is true
     *  (the default), write a warning message to the console when the requested
     *  property is not available.
     */
    double get(const std::string& propertyName, const std::string& key, double alt, bool warn=true) const;
    /** Set property \a propertyName of item \a key to value. If that property
     *  already exists, its value will be modified, otherwise the property will
     *  be added.
     */
    void set(const std::string& propertyName, const std::string& key, double value);
    /** Return a constant reference to the json object in which all properties
     *  are stored.
     */
    const json_type& data() const { return m_data; }
private:
    /** Read properties from file \a fname. That file must be a LoKI-B
     *  database file. Percentage signs start comments, and lines that are
     *  non-empty after stripping the comments must contain key-value pairs,
     *  the values must be valid double values.
     *  The properties are stored as elements of a JSON object that contains
     *  the values for all keys and is returned by this function. A runtime
     *  error is thrown when the file cannot be read or contains bad data.
     */
    json_type readGasPropertyFile(const std::filesystem::path& fname) const;
    /// The json object in which all properties are stored.
    json_type m_data;
};

} // namespace loki

#endif // LOKI_CPP_GASPROPERTIES_H
