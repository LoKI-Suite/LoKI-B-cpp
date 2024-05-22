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

void patchStateProperty(json_type& stateProp)
{
    // make a copy, clear fracs and re-add the properties in new form
    const std::vector<std::string> entryVector = stateProp;
    stateProp = json_type{};
    for (const auto &line : entryVector)
    {
        // look for a line of the form "S = E"
        static const std::regex reProperty(R"((\S+)\s+=\s+(\S+)\s*$)");
        std::smatch m;
        if (std::regex_search(line, m, reProperty))
        {
            // Found. E can be a literal value, or a function with arguments.
            // value or function
            const std::string state_id = m.str(1);
            const std::string expr = m.str(2);

            // 2. Now apply the expression.
            // Try to parse expr as a number first...
            double value;
            if (Parse::getValue(expr, value))
            {
                json_type prop;
                prop["type"] = "constant";
                prop["value"] = value;
                stateProp["states"][state_id] = prop;
            }
            else
            {
                // expr is not a number. We treat is a function (maybe with arguments).
                static const std::regex reFuncArgs(R"(\s*(\w+)@?(.*))");
                std::smatch fm;
                if (!std::regex_match(expr, fm, reFuncArgs))
                {
                    throw std::runtime_error("Could not parse function "
                        "name and argument list from string '"
                        + expr + "'.");
                }
                const std::string functionName = fm.str(1);
                const std::string argumentString = fm.str(2);

                json_type prop;
                prop["type"] = "function";
                prop["name"] = functionName;

                // create an argument list for the function (possibly empty)
                std::vector<double> arguments;
                static const std::regex reArgList(R"(\s*([\w\.]+)\s*(?:[,\]]|$))");
                for (auto it = std::sregex_iterator(argumentString.begin(), argumentString.end(), reArgList);
                     it != std::sregex_iterator(); ++it)
                {
                    // pvalue set by getValue
                    json_type val_node;
                    const std::string arg{it->str(1)};
                    double pvalue;
                    if (Parse::getValue(arg, pvalue))
                    {
                        val_node = pvalue;
                    }
                    else
                    {
                        val_node = arg;
                    }
                    prop["arguments"].push_back(val_node);
                }
                stateProp["states"][state_id] = prop;
            }
        }
        else
        {
            stateProp["files"].push_back(line);
        }
    }
}

void patchStateProperties(json_type& props)
{
    for (auto& prop : props)
    {
        patchStateProperty(prop);
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
    patchStateProperties(ek.at("stateProperties"));
    if (ek.contains("effectiveCrossSectionPopulations"))
    {
        json_type effPop;
        for (const auto& fileName : ek.at("effectiveCrossSectionPopulations"))
        {
            effPop["files"].push_back(fileName.get<std::string>());
        }
        ek.at("effectiveCrossSectionPopulations") = effPop;
    }
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

json_type readLegacyGasPropertyFile(const std::filesystem::path& fname)
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

json_type readLegacyStatePropertyFile(const std::filesystem::path &fileName)
{
    Log<Message>::Notify("Processing state property file '" + fileName.generic_string() + "'.");
    json_type entries;
    /** \bug 'S 1.2.3' will be accepted by this regex. Subsequently,
     *        getValue will result in the value 1.2, since that does
     *        not care about trailing characters.
     */
    static const std::regex word(R"((\S+)\s*$)");
    static const std::regex assign(R"((\S+)\s+(\S+)\s*$)");

    std::string fileBuffer;
    if (!Parse::stringBufferFromFile(fileName, fileBuffer))
    {
        Log<Message>::Error("Could not open state property file '" + fileName.generic_string() + "' for reading.");
    }
    std::stringstream ss{fileBuffer};
    std::string line;
    while (std::getline(ss, line))
    {
        std::smatch m;
        if (std::regex_search(line, m, assign))
        {
            const std::string stateString = m.str(1);
            const std::string valueString = m.str(2);
            entries.push_back(stateString + " = " + valueString);
        }
        else if (std::regex_search(line, m, word))
        {
            const std::string fname = m.str(1);
            entries.push_back(fname);
        }
        else
        {
            throw std::runtime_error("Syntax error in file '" + fileName.generic_string() + "', line '" + line +
                                     "': expected '<states> = <value|function call>' or a file name.");
        }
    }
    patchStateProperty(entries);
    return entries;
}

} // namespace loki
