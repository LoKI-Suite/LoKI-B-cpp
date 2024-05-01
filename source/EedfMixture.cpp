//
// Created by daan on 21-5-19.
//

#include "LoKI-B/EedfMixture.h"
#include "LoKI-B/Log.h"

#include <stdexcept>

namespace loki
{

EedfMixture::EedfMixture(const std::filesystem::path &basePath, const Grid *grid, const json_type &cnf, const WorkingConditions *workingConditions)
{
    if (cnf.contains("mixture"))
    {
        m_collision_data.loadCollisionsJSON(cnf.at("mixture"), m_composition, grid, false);
    }
    if (cnf.contains("LXCatFiles"))
    {
        loadCollisions(basePath, cnf.at("LXCatFiles").get<std::vector<std::string>>(), grid);
    }
    if (cnf.contains("LXCatFilesExtra"))
    {
        loadCollisions(basePath, cnf.at("LXCatFilesExtra").get<std::vector<std::string>>(), grid, true);
    }

    this->loadGasProperties(basePath, cnf.at("gasProperties"));
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
    composition().loadStateProperties(basePath, cnf.at("stateProperties"), workingConditions);
    composition().evaluateReducedDensities();

    for (auto &cd : m_collision_data.data_per_gas())
    {
        if (!cd.isDummy())
        {
            cd.checkElasticCollisions(electron, grid);
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

void EedfMixture::loadCollisions(const std::filesystem::path &basePath, const std::vector<std::string> &files, const Grid *energyGrid, bool isExtra)
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
            m_collision_data.loadCollisionsJSON(cnf, m_composition, energyGrid, isExtra);
        }
        else
        {
            m_collision_data.loadCollisionsClassic(fileName, m_composition, energyGrid, isExtra);
        }
    }
    Log<Message>::Notify("Finished loading collisions.");
}

void EedfMixture::loadGasProperties(const std::filesystem::path &basePath, const json_type &cnf)
{
    composition().loadGasProperties(basePath, cnf);
    readGasProperty(basePath, composition().gases(), cnf, "OPBParameter", false,
        [this](Gas& gas, double value) { this->collision_data().get_data_for(gas).setOPBParameter(value); } );
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
