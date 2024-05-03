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

#include "LoKI-B/GasProperties.h"
#include "LoKI-B/Constant.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/Parse.h"

#include <string>

namespace loki {

GasProperties::GasProperties()
{
    // If not yet available, add the "mass" property of "e", the electron.
    if (!has("mass","e"))
    {
        set("mass","e",Constant::electronMass);
    }
}

GasProperties::GasProperties(const std::filesystem::path &basePath, const json_type& pnode)
{
    for (const auto& p : pnode.items())
    {
        if (p.value().type() == json_type::value_t::string)
        {
            // We have, for example, { "mass": "DataBases/masses.txt" }
            // Replace this with a member "mass: " [ { "gasname", "mass" }, ... ]
            std::filesystem::path path(p.value());
            if (path.is_relative())
            {
                path = basePath.parent_path() / path;
            }
            m_data[p.key()] = readGasPropertyFile(path);
        }
        else if (p.value().type() == json_type::value_t::object)
        {
            // Assume that the data are already in the correct form.
            // Just copy them in the correct place.
            m_data[p.key()] = p.value();
        }
        else
        {
            std::cout << "Gas properties: skipping: " << p.key() << ": " << p.value().dump(2) << std::endl;
        }
    }
    // If not yet available, add the "mass" property of "e", the electron.
    if (!has("mass","e"))
    {
        set("mass","e",Constant::electronMass);
    }
}

json_type GasProperties::readGasPropertyFile(const std::filesystem::path& fname) const
{
    if (fname.extension()==".json")
    {
        return read_json_from_file(fname);
    }
    json_type result;
    std::string str;
    if (!Parse::stringBufferFromFile(fname,str))
    {
        throw std::runtime_error("Error opening/reading file '"
            + fname.generic_string() + "'.");
    }
    std::stringstream ss;
    ss << str;
    while (!ss.eof())
    {
        std::string gas;
        double value;
        ss >> gas >> value;
        if (ss.fail())
        {
            throw std::runtime_error("Bad data in file '"
                + fname.generic_string() + "'.");
        }
        result[gas] = value;
    }
    return result;
}

bool GasProperties::has(const std::string& propertyName) const
{
    return data().contains(propertyName);
}

bool GasProperties::has(const std::string& propertyName, const std::string& key) const
{
    return data().contains(propertyName) && data()[propertyName].contains(key);
}

double GasProperties::get(const std::string& propertyName, const std::string& key) const
{
    return data().at(propertyName).at(key);
}

double GasProperties::get(const std::string& propertyName, const std::string& key, double alt, bool warn) const
{
    const json_type* res_ptr = data().contains(propertyName)
        ? (data()[propertyName].contains(key) ? &data()[propertyName][key] : nullptr)
        : nullptr;
    if (res_ptr)
    {
        return *res_ptr;
    }
    if (warn)
    {
        Log<GasPropertyError>::Warning(propertyName + " for " + key);
    }
    return alt;
}

void GasProperties::set(const std::string& propertyName, const std::string& key, double value)
{
    m_data[propertyName][key] = value;
}

} // namespace loki
