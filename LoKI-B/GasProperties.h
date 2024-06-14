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
 *  \date   2 May 2024
 */

#ifndef LOKI_CPP_GASPROPERTIES_H
#define LOKI_CPP_GASPROPERTIES_H

#include "LoKI-B/json.h"
#include <filesystem>

namespace loki {

/** This class manages the information that is contained in a LoKI-B
 *  'gasProperties' section. It stores the data in a json object. The class
 *  provides members \c has() to check if particular data are available,
 *  \c get() to retrieve property values and \c set to add or modify properties.
 *
 *  Member \c data() returns a constant reference to the entire data collection.
 *  The object is structured as a map-of-maps. The keys of the elements are
 *  the names of the properties that this class knows about, the values are
 *  objects themselves, consisting of pairs of (gas) names and property values.
 *  This is illustrated by the following example:
 *
 *  \verbatim
    {
        "mass": {
            "N2": 4.651834066656e-26
            "O2": 5.313392820191999e-26,
            "e": 9.10938356e-31,
            ...
        },
        "rotationalConstant": {
            "N2": 2.477501826053967e-4,
            "O2": 1.782496009128667e-4
            ...
        },
        ...
    } \endverbatim
 *
 *  \author Jan van Dijk
 *  \date   2 May 2024
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
     *  found, a runtime_error is thrown.
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
    /// The json object in which all properties are stored.
    json_type m_data;
};

} // namespace loki

#endif // LOKI_CPP_GASPROPERTIES_H
