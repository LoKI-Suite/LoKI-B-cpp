//
// Created by daan on 21-5-19.
//

#include "LoKI-B/EedfMixture.h"
#include "LoKI-B/Log.h"

#include <stdexcept>

namespace loki
{

EedfMixture::EedfMixture(Grid *grid, const ElectronKineticsSetup &setup,
                               const WorkingConditions *workingConditions)
{
    loadCollisions(setup.LXCatFiles, grid);
    if (!setup.LXCatFilesExtra.empty())
    {
        loadCollisions(setup.LXCatFilesExtra, grid, true);
    }

    this->loadGasProperties(setup.gasProperties);
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
    /** For the electron, the root state has a population of 0 (the default),
     *  other states in the hierarchy (which have no siblings) 1. This code
     *  can be simplified once a final decision has been made how to represent
     *  the electron State. As a charged state? An electronic state? Or just as
     *  the root of the gas? At present, a charged state seems best, because the
     *  charge of the electron gas is meaningful. If we move to a situation where
     *  ionic varieties (such as N2^+) are made separate gases, and the gas itself
     *  hosts the properties (such as charge), it makes more sens to let the
     *  root be 'the' electronic state, since there is no need for additional structure.
     */
    // for (State *s = electron; s->type() != StateType::root; s = s->parent())
    // {
    //     s->setPopulation(1.0);
    // }
    composition().loadStateProperties(setup.stateProperties, workingConditions);
    composition().evaluateReducedDensities();

    for (auto& cd : m_collision_data.data_per_gas())
    {
        if (!cd.isDummy())
        {
            cd.checkElasticCollisions(electron, grid);
        }
    }

    for (const std::string& gasName : setup.CARgases)
    {
        this->addCARGas(gasName);
    }

    //        this->evaluateTotalAndElasticCS();
}

EedfMixture::EedfMixture(Grid *grid, const json_type &cnf, const WorkingConditions *workingConditions)
{
    if (cnf.contains("mixture"))
    {
        m_collision_data.loadCollisionsJSON(cnf.at("mixture"), m_composition, grid, false);
    }
    if (cnf.contains("LXCatFiles"))
    {
        loadCollisions(cnf.at("LXCatFiles").get<std::vector<std::string>>(), grid);
    }
    if (cnf.contains("LXCatFilesExtra"))
    {
        loadCollisions(cnf.at("LXCatFilesExtra").get<std::vector<std::string>>(), grid, true);
    }

    this->loadGasProperties(cnf.at("gasProperties"));
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
    composition().loadStateProperties(cnf.at("stateProperties"), workingConditions);
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

void EedfMixture::loadCollisions(const std::vector<std::string> &files, Grid *energyGrid, bool isExtra)
{
    Log<Message>::Notify("Started loading collisions.");
#ifdef EMSCRIPTEN
    const std::string inputPath{""};
#else
    const std::string inputPath{"../Input/"};
#endif

    for (const std::string &f : files)
    {
        const std::string fname = inputPath + f;
        if (fname.size() >= 5 && fname.substr(fname.size() - 5) == ".json")
        {
            const json_type cnf = read_json_from_file(fname);
            m_collision_data.loadCollisionsJSON(cnf, m_composition, energyGrid, isExtra);
        }
        else
        {
            m_collision_data.loadCollisionsClassic(fname, m_composition, energyGrid, isExtra);
        }
    }
    Log<Message>::Notify("Finished loading collisions.");
}

void EedfMixture::loadGasProperties(const GasPropertiesSetup &setup)
{
    composition().loadGasProperties(setup);
    readGasPropertyFile(composition().gases(), setup.OPBParameter, "OPBParameter", false,
        [this](Gas& gas, double value) { this->collision_data().get_data_for(gas).setOPBParameter(value); } );
}

void EedfMixture::loadGasProperties(const json_type &cnf)
{
    composition().loadGasProperties(cnf);
    readGasProperty(composition().gases(), cnf, "OPBParameter", false,
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
