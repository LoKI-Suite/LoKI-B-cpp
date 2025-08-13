/** \file
 *
 *  Implementation of a class that managed the properties of a mixture.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2025 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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
 *  \date   21 May 2019
 */

#include "LoKI-B/EedfMixture.h"
#include "LoKI-B/GasProperties.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/LegacyToJSON.h"

#include <stdexcept>

namespace loki
{

EedfMixture::EedfMixture(const std::filesystem::path &basePath, const Grid *grid, const json_type &cnf, const WorkingConditions *workingConditions)
{
    const GasProperties gasProps(basePath, cnf.at("gasProperties"));
    if (cnf.contains("mixture"))
    {
        m_collision_data.loadCollisionsJSON(cnf.at("mixture"), gasProps, m_composition, grid, false);
    }
    if (cnf.contains("LXCatFiles"))
    {
        loadCollisions(basePath, cnf.at("LXCatFiles").get<std::vector<std::string>>(), gasProps, grid);
    }
    if (cnf.contains("LXCatFilesExtra"))
    {
        loadCollisions(basePath, cnf.at("LXCatFilesExtra").get<std::vector<std::string>>(), gasProps, grid, true);
    }

    for (auto &gas : composition().gases())
    {
        gas->propagateFraction();
    }

    /// \todo Store the electron state pointer as a member? Does the state type matter?
    Gas::State *electron = m_composition.findState(StateEntry::electronEntry());
    if (!electron)
    {
        throw std::runtime_error("No electron species has been found in the process list.");
    }
    /// \todo See the comments about the population in the other overload
    m_composition.loadStateProperties(basePath, cnf.at("stateProperties"), workingConditions);
    m_composition.evaluateReducedDensities();

    EffectivePopulationsMap effectivePopulations;
    if (cnf.contains("effectiveCrossSectionPopulations"))
    {
        readEffectivePopulations(basePath,cnf.at("effectiveCrossSectionPopulations"),effectivePopulations);
    }
    for (auto &cd : m_collision_data.data_per_gas())
    {
        if (!cd.isDummy())
        {
            cd.checkElasticCollisions(electron, grid, effectivePopulations);
        }
    }
    if (cnf.contains("CARgases"))
    {
        for (const auto& cg : cnf.at("CARgases"))
        {
            addCARGas(cg.get<std::string>());
        }
    }

    //        this->evaluateTotalAndElasticCS();
}

void EedfMixture::readEffectivePopulations(const std::filesystem::path &basePath, const json_type& effPop, EffectivePopulationsMap& effectivePopulations) const
{
    if (effPop.contains("states"))
    {
        for (const auto& e : effPop.at("states").items())
        {
            const StateEntry entry = propertyStateFromString(e.key());
            Gas::State::ChildContainer states = const_cast<GasMixture&>(composition()).findStates(entry);

            for (const auto &state : states) {
                if (effectivePopulations.count(state))
                {
                    throw std::runtime_error("Duplicate specification for state '"
                        + e.key() + "' found.");
                }
                if (e.value().at("type")=="constant")
                {
                    effectivePopulations[state] = e.value().at("value");
                    std::cout << "readEffectivePopulations: setting effective cross section "
                                 "population of state '" << e.key() << " to "
                              << e.value().at("value").dump(2) << "'." << std::endl;
                }
                else
                {
                    /// \todo Should we also support functions here?
                    throw std::runtime_error("readEffectivePopulations for state '"
                        + e.key() + "': only type 'Constant' is supported.");
                }
            }
        }
    }
    if (effPop.contains("files"))
    {
        for (const auto& f : effPop.at("files"))
        {
            try {
                std::cout << "readEffectivePopulations: handling file '" + f.get<std::string>() + "'." << std::endl;
                readEffectivePopulations(basePath,f.get<std::string>(),effectivePopulations);
            }
            catch(std::exception& exc)
            {
                throw std::runtime_error("While reading '" + f.get<std::string>() + ": " + exc.what());
            }
        }
    }
}

void EedfMixture::readEffectivePopulations(const std::filesystem::path &basePath, const std::string& f, EffectivePopulationsMap& effectivePopulations) const
{
    std::filesystem::path fileName(f);

    if (fileName.is_relative()) {
        fileName = basePath.parent_path() / fileName;
    }

    json_type effPop;
    if (fileName.has_extension() && fileName.extension() == ".json")
    {
        effPop = read_json_from_file(fileName);
    }
    else
    {
        effPop = readLegacyStatePropertyFile(fileName);
    }
    readEffectivePopulations(basePath,effPop,effectivePopulations);
}

void EedfMixture::loadCollisions(const std::filesystem::path &basePath, const std::vector<std::string> &files, const GasProperties& gasProps, const Grid *energyGrid, bool isExtra)
{
    Log<Message>::Notify("Started loading collisions.");

    for (const std::string &f : files)
    {
        std::filesystem::path fileName(f);

        if (fileName.is_relative()) {
            fileName = basePath.parent_path() / fileName;
        }

        if (fileName.has_extension() && fileName.extension() == ".json")
        {
            const json_type cnf = read_json_from_file(fileName);
            m_collision_data.loadCollisionsJSON(cnf, gasProps, m_composition, energyGrid, isExtra);
        }
        else
        {
            m_collision_data.loadCollisionsClassic(fileName, gasProps, m_composition, energyGrid, isExtra);
        }
    }
    Log<Message>::Notify("Finished loading collisions.");
}

void EedfMixture::addCARGas(const std::string& gasName)
try {
    const Gas *gas = composition().findGas(gasName);
    if (gas == nullptr)
    {
        throw std::runtime_error("Gas not found.");
    }
    if (gas->electricQuadrupoleMoment < 0)
    {
        throw std::runtime_error("Electrical quadrupole moment not specified.");
    }
    if (gas->rotationalConstant < 0)
    {
        throw std::runtime_error("Rotational constant not specified.");
    }
    if (!m_collision_data.get_data_for(*gas).collisions(CollisionType::rotational).empty())
    {
        throw std::runtime_error("Rotational collision have been specified.");
    }
    // all checks passed. add the gas to the list of CAR gases
    m_CARGases.emplace_back(gas);
    Log<Message>::Notify("Registered gas '" + gasName + "' as CAR gas.");
}
catch(std::exception& exc)
{
    throw std::runtime_error("Error adding CAR gas '" + gasName + "': "
        + std::string{exc.what()});
}

} // namespace loki
