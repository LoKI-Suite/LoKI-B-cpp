/** \file
 *
 *  Implementation of the GasMixture class.
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
 *  \author Daan Boer and Jan van Dijk (C++ version)
 *  \date   May 2019
 */

#include "LoKI-B/GasMixture.h"
#include "LoKI-B/LegacyToJSON.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/PropertyFunctions.h"

#include <filesystem>

namespace loki
{

GasMixture::~GasMixture()
{
}

Gas *GasMixture::addGas(const GasProperties& gasProps, const std::string& name)
{
    if (findGas(name))
    {
        throw std::logic_error("Attempt to register gas with name '" + name + "' twice.");
    }
    return m_gases.emplace_back(new Gas(gasProps,name)).get();
}

Gas *GasMixture::ensureGas(const GasProperties& gasProps, const std::string &name)
{
    Gas *gas = findGas(name);
    if (!gas)
    {
        gas = addGas(gasProps,name);
    }
    return gas;
}

Gas::State *GasMixture::ensureState(const GasProperties& gasProps, const StateEntry &entry)
{
    typename StateMap::iterator it = m_states.find(entry.m_id);
    if (it != m_states.end())
    {
        return it->second;
    }
    Gas::State *state = ensureGas(gasProps,entry.m_gasName)->ensureState(entry);
    m_states[entry.m_id] = state;
    return state;
}

Gas::State *GasMixture::findStateById(const std::string &stateId)
{
    typename StateMap::iterator it = m_states.find(stateId);
    return it == m_states.end() ? nullptr : it->second;
}
void GasMixture::print(std::ostream &os)
{
    for (const auto &gas : m_gases)
    {
        os << "Gas: " << gas->name() << std::endl;
        gas->print(os);
    }
}

void GasMixture::checkGasFractions() const
{
    double norm = 0;

    for (const auto &gas : m_gases)
    {
        norm += gas->fraction;
    }

    if (std::abs(norm - 1.) > 10. * std::numeric_limits<double>::epsilon())
        Log<Message>::Error("Gas fractions are not properly normalized.");
}

void GasMixture::checkPopulations() const
{
    for (auto &gas : m_gases)
        gas->checkPopulations();
}

const Gas *GasMixture::findGas(const std::string &name) const
{
    auto it = std::find_if(m_gases.begin(), m_gases.end(),
                           [&name](const std::unique_ptr<Gas> &gas) { return gas->name() == name; });
    return it == m_gases.end() ? nullptr : it->get();
}

Gas *GasMixture::findGas(const std::string &name)
{
    auto it = std::find_if(m_gases.begin(), m_gases.end(),
                           [&name](const std::unique_ptr<Gas> &gas) { return gas->name() == name; });
    return it == m_gases.end() ? nullptr : it->get();
}

Gas::State *GasMixture::findState(const StateEntry &entry)
{
    Gas *gas = findGas(entry.m_gasName);
    return gas ? gas->findState(entry) : nullptr;
}

Gas::State::ChildContainer GasMixture::findStates(const StateEntry &entry)
{
    using ChildContainer = Gas::State::ChildContainer;
    Gas::State* state = findState(entry);
    return (!state) ? ChildContainer{} : entry.hasWildCard() ? state->siblings() : ChildContainer{ state };
}

void GasMixture::loadStatePropertyEntry(const std::string& state_id, const json_type& propEntry,
                                   StatePropertyType propertyType, const WorkingConditions *workingConditions)
{
    std::cout << "loadStatePropertyEntry: handling states " << state_id << ":\n" << propEntry.dump(2) << std::endl;
    // 1. Find the state(s) for which the property must be set.
    const StateEntry entry = propertyStateFromString(state_id);
    /** \todo Should this be part of propertyStateFromString?
     *        Is there a reason to accept a 'none'-result?
     */
    if (entry.m_level == none)
    {
        throw std::runtime_error("loadStateProperty: illegal "
                        "state identifier '" + propEntry.at("states").get<std::string>() + "'.");
    }
    Gas::State::ChildContainer states = findStates(entry);
    if (states.empty())
    {
        throw std::runtime_error("loadStateProperty: could not find "
                        "state or state group '" + propEntry.at("states").get<std::string>() + "'.");
    }

    // 2. Now apply the expression.
    if (propEntry.at("type")=="constant")
    {
        PropertyFunctions::constantValue(states, propEntry["value"], propertyType);
    }
    else if (propEntry.at("type")=="function")
    {
        const std::string functionName = propEntry.at("name");
        // create an argument list for the function (possibly empty)
        std::vector<double> arguments;
        if (propEntry.contains("arguments"))
        {
            for (const auto& arg : propEntry.at("arguments"))
            {
                double pvalue;
                if (arg.is_number())
                {
                    pvalue = arg;
                }
                else if (arg.is_string())
                {
                    pvalue = workingConditions->getParameter(arg);
                }
                else
                {
                    throw std::runtime_error("Argument '" + arg.dump(2) + "' is neither a numerical "
                                "value, nor a parameter name.");
                }
                arguments.emplace_back(pvalue);
            }
        }
        PropertyFunctions::callByName(functionName, states, arguments, propertyType);
    }
    else
    {
        throw std::runtime_error("Unknown property type '" + propEntry.at("type").get<std::string>() + "'.");
    }
}

void GasMixture::loadStateProperty(const std::filesystem::path &basePath, const json_type& stateProp,
                                   StatePropertyType propertyType, const WorkingConditions *workingConditions)
{
    if (stateProp.contains("states"))
    {
        for (const auto& propEntry : stateProp["states"].items())
        {
            loadStatePropertyEntry(propEntry.key(),propEntry.value(), propertyType, workingConditions);
        }
    }
    if (stateProp.contains("files"))
    {
        for (const auto &propEntry : stateProp["files"])
        {
            std::cout << "loadStateProperty: handling: " << propEntry.dump() << std::endl;
            std::filesystem::path fileName(propEntry);
            if (fileName.is_relative())
            {
                fileName = basePath.parent_path() / fileName;
            }
            try
            {
                // call this function recursively for each specified file.
                if (fileName.extension()==".json")
                {
                    const json_type entries = read_json_from_file(fileName);
                    loadStateProperty(basePath, entries, propertyType, workingConditions);
                }
                else
                {
                    // support for legacy (".in") files.
                    const json_type entries = readLegacyStatePropertyFile(fileName);
                    loadStateProperty(basePath, entries, propertyType, workingConditions);
                }
            }
            catch (std::exception& exc)
            {
                throw std::runtime_error("Error configuring state property '"
                                        + statePropertyName(propertyType) + "': "
                                        + std::string{exc.what()} );
            }
        }
    }
}

void GasMixture::evaluateReducedDensities()
{
    for (auto &gas : m_gases)
    {
        gas->evaluateReducedDensities();
    }
}

void GasMixture::loadStateProperties(const std::filesystem::path &basePath, const json_type &cnf, const WorkingConditions *workingConditions)
{
    if (cnf.contains("energy"))
    {
        loadStateProperty(basePath, cnf.at("energy"), StatePropertyType::energy, workingConditions);
    }
    if (cnf.contains("statisticalWeight"))
    {
        loadStateProperty(basePath, cnf.at("statisticalWeight"), StatePropertyType::statisticalWeight, workingConditions);
    }
    // the population section is mandatory
    loadStateProperty(basePath, cnf.at("population"), StatePropertyType::population, workingConditions);

    checkPopulations();
}

} // namespace loki
