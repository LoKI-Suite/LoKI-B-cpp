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

#include "LoKI-B/Log.h"
#include "LoKI-B/OffSideToJSON.h"
#include "LoKI-B/Parse.h"
#include <regex>

namespace loki
{

void patchRange(json_type& cnf, const std::string& unit)
{
	//std::cout << "patchRange, unit = " << unit << ". Before: " << cnf.dump(2) << std::endl;
	if (cnf.is_string())
	{
		cnf = Parse::makeJsonFromFunctionCall(cnf.get<std::string>(),
			"range", {"min","max","nvalues"});
	}
	else if (cnf.is_number())
	{
		json_type cnf2;
		cnf2["value"] = cnf;
		cnf2["unit"] = unit;
		cnf = cnf2;
	}
	else
	{
		/// \todo Check that this is all we need to check for
		// nothing
	}
	//std::cout << "patchRange, after: " << cnf.dump(2) << std::endl;
}

void patchWorkingConditions(json_type& wc)
{
	const std::vector<std::pair<std::string,std::string>> params
	{
		{"gasPressure", "Pa"},
		{"gasTemperature", "K"},
		{"electronDensity", "m^-3"},
		{"electronTemperature", "eV"},
		{"excitationFrequency", "Hz"},
		{"reducedField", "Td"}
	};
	for (const auto& p : params)
	{
		if (wc.contains(p.first))
		{
			patchRange(wc[p.first],p.second);
		}
	}	
}

/* change an array of strings of the form "X = 1.0" into a object with
 * elements of the form '"X": 1.0'
 */
void patchFractions(json_type& fracs)
{
    json_type new_fracs;
    // Parse fractions
    const std::regex r(R"(([\w\d]*)\s*=\s*(\d*\.?\d*)$)");
    std::smatch m;

    for (json_type& frac : fracs)
    {
        try
        {
            const std::string propertyString = frac;
            if (!std::regex_search(propertyString, m, r))
            {
                Log<Message>::Error("Could not parse gas fraction '" + frac.get<std::string>() + "'.");
            }
            const auto name = m.str(1);
            const double value = Parse::getValue(m.str(2));
            new_fracs[m.str(1)] = value;
        }
        catch (std::exception& exc)
        {
            Log<Message>::Error("Could not parse gas fraction '" + frac.get<std::string>() + "': " + exc.what());
        }
    }
    fracs = new_fracs;
}


void patchElectronKinetics(json_type& ek)
{
	patchFractions(ek.at("gasProperties").at("fraction"));
}

json_type legacyToJSON(const json_type& legacy)
{
	json_type json(legacy);
	if (json.contains("workingConditions"))
	{
		patchWorkingConditions(json["workingConditions"]);
	}
	if (json.contains("electronKinetics"))
	{
		patchElectronKinetics(json["electronKinetics"]);
	}
	return json;
}

json_type legacyToJSON(std::istream& is)
{
	return legacyToJSON(offSideToJSON(is));
}

json_type legacyToJSON(const std::filesystem::path& fname)
{
	return legacyToJSON(offSideToJSON(fname));
}

} // namespace loki
