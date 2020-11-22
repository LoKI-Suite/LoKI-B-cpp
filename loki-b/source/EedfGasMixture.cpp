//
// Created by daan on 20-5-19.
//

#include "LoKI-B/EedfGasMixture.h"
#include "LoKI-B/Constant.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/json.h"
#include <iostream>
#include <stdexcept>
#include <regex>

namespace loki
{

EedfGasMixture::EedfGasMixture(Grid *grid, const ElectronKineticsSetup &setup,
                               const WorkingConditions *workingConditions)
    : grid(grid)
{

    loadCollisions(setup.LXCatFiles, grid);

    if (!setup.LXCatFilesExtra.empty())
        loadCollisions(setup.LXCatFilesExtra, grid, true);

    this->loadGasProperties(setup.gasProperties);

    for (auto &gas : gases())
        gas->propagateFraction();

    /// \todo Store the electron state pointer as a member? Does the state type matter?
    GasBase::State *electron = ensureState(StateEntry::electronEntry());
    /** For the electron, the root state has a population of 0 (the default),
     *  other states in the hierarchy (which have no siblings) 1. This code
     *  can be simplified once a final decision has been made how to represent
     *  the electron State. As a charged state? An electronic state? Or just as
     *  the root of the gas? At present, a charged state seems best, because the
     *  charge of the electron gas is meaninful. If we move to a situation where
     *  ionic varieties (such as N2^+) are made separate gases, and the gas itself
     *  hosts the properties (such as charge), it makes more sens to let the
     *  root be 'the' electronic state, since there is no need for additional structure.
     */
    // for (State *s = electron; s->type != StateType::root; s = s->parent())
    // {
    //     s->population = 1.0;
    // }
    this->loadStateProperties(setup.stateProperties, workingConditions);
    this->evaluateStateDensities();

    for (auto &gas : gases())
    {
        if (!gas->isDummy())
            gas->checkElasticCollisions(electron, grid);
    }

    this->addCARGases(setup.CARgases);

    //        this->evaluateTotalAndElasticCS();
}

EedfGasMixture::EedfGasMixture(Grid *grid, const json_type &cnf, const WorkingConditions *workingConditions)
    : grid(grid)
{

    loadCollisions(cnf.at("LXCatFiles"), grid);

    if (cnf.contains("LXCatFilesExtra"))
        loadCollisions(cnf.at("LXCatFilesExtra").get<std::vector<std::string>>(), grid, true);

    this->loadGasProperties(cnf.at("gasProperties"));

    for (auto &gas : gases())
        gas->propagateFraction();

    /// \todo Store the electron state pointer as a member? Does the state type matter?
    GasBase::State *electron = ensureState(StateEntry::electronEntry());
    /**
     */
    // for (State *s = electron; s->type != StateType::root; s = s->parent())
    // {
    //     s->population = 1.0;
    // }
    this->loadStateProperties(cnf.at("stateProperties"), workingConditions);
    this->evaluateStateDensities();

    for (auto &gas : gases())
    {
        if (!gas->isDummy())
            gas->checkElasticCollisions(electron, grid);
    }
    if (cnf.contains("CARgases"))
    {
        this->addCARGases(cnf.at("CARgases"));
    }

    //        this->evaluateTotalAndElasticCS();
}

void EedfGasMixture::loadCollisions(const std::string &file, Grid *energyGrid, bool isExtra)
{
    const std::regex reParam(R"(PARAM\.:)");
    /** \todo Valid doubles like 1.57e1 are not matched, in which case the threshold
     *        specification will be ignored and 0 will be assumed.
     */
    const std::regex reThreshold(R"(E = (\d*\.?\d*) eV)");
    const std::regex reProcess(R"(\[(.+?)(<->|->)(.+?), (\w+)\])");
    std::ifstream in(file);
    if (!in.is_open())
    {
        Log<FileError>::Warning(file);
        return;
    }

    std::string line;
    std::smatch mThreshold, mProcess;

    /*  Expectations: we cycle throught the file and look for a line that has "PARAM.:".
     *  When we find that, we start parsing a process description.
     *  1. Read the arguments of the parameter line, see if a threshold is specified ("E = <double> eV"),
     *     Set the threshold to 0 otherwise.
     *  2. The next line should describe the process: it must be of the form [<lhs> -> <rhs>, <type>].
     *  3. A Collision is created from the process description information
     *  4. (Maybe) a cross section object is attached to the collision object, based on subsequent lines.
     *
     */
    while (std::getline(in, line))
    {
        try
        {
            if (!std::regex_search(line, reParam))
            {
                continue;
            }
            double threshold = 0.0;
            if (std::regex_search(line, mThreshold, reThreshold))
            {
                std::stringstream ss(mThreshold[1]);
                if (!(ss >> threshold))
                    Log<Message>::Error("Non numerical threshold value.");
            }
            if (!std::getline(in, line) || !std::regex_search(line, mProcess, reProcess))
            {
                Log<LXCatError>::Error(file);
            }

            try
            {
                // mProcess[1...4] represent: lhs, separator, rhs and type, and could
                // be something like:
                //  1: "He + e"
                //  2: "->" ("<->" for a two-way process)
                //  3: "He + e"
                //  4: "Elastic"

                std::vector<StateEntry> entry_lhsStates, entry_rhsStates;
                std::vector<uint16_t> entry_lhsCoeffs, entry_rhsCoeffs;

                entriesFromString(mProcess[1].str(), entry_lhsStates, &entry_lhsCoeffs);
                const bool reverseAlso = (mProcess[2].str()[0] == '<');
                entriesFromString(mProcess[3].str(), entry_rhsStates, &entry_rhsCoeffs);
                const CollisionType entry_type = getCollisionType(mProcess[4].str());

                // 1. Create vectors of pointers to the states that appear
                //    on the left and right-hand sides of the process.
                //    Create the states when necessary.
                std::vector<State *> lhsStates;
                std::vector<State *> rhsStates;
                for (auto &stateEntry : entry_lhsStates)
                {
                    lhsStates.emplace_back(ensureState(stateEntry));
                }
                for (auto &stateEntry : entry_rhsStates)
                {
                    rhsStates.emplace_back(ensureState(stateEntry));
                }
                // 2. Create the collision object (but do not configure a CrossSection
                //    object yet).
                std::unique_ptr<Collision> collision{
                    new Collision(entry_type, lhsStates, entry_lhsCoeffs, rhsStates, entry_rhsCoeffs, reverseAlso)};
                Log<Message>::Notify(*collision);
                // 3. See if we already have a collision of the same type with the same lhs and rhs.
                //    If we do, isue a warning, discard the newly created collision object and return
                //    nullptr. (NOTE: collision is a unique_ptr and will delete the collision object,
                //    since we do not release that before returning from this function..)
                for (const auto &c : m_collisions)
                {
                    if (c->is_same_as(*collision))
                    {
                        Log<DoubleCollision>::Warning(*c);
                        collision.release();
                        break;
                    }
                }
                // 4. Register the collision with the (single) gas that is the target of the
                //    binary electron-heavy process, add it to out own list of collisions
                //    and return the (released) pointer to our caller.
                //    That the left-hand side of the process is indeed of the form 'e + X' is
                //    checked by the EedfCollision constructor, and the target 'X' is returned
                //    by its member getTarget().
                /** \todo Eliminate the following cast at some moment: we can do
                 *        EedfGas& eedfGas = *ensureGas(collision->getTarget()->gas().name);
                 *        but could check/assert that the pointer is the same.
                 */
                if (collision.get())
                {
                    EedfGas &eedfGas = static_cast<EedfGas &>(collision->getTarget()->gas());
                    eedfGas.addCollision(collision.get(), isExtra);
                    m_collisions.push_back(collision.get());
                    hasCollisions[static_cast<uint8_t>(collision->type)] = true;

                    const bool isElasticOrEffective =
                        (collision->type == CollisionType::effective || collision->type == CollisionType::elastic);
                    collision->crossSection.reset(new CrossSection(threshold, energyGrid, isElasticOrEffective, in));
                    collision.release();
                }
            }
            catch (std::exception &exc)
            {
                throw std::runtime_error("Error while parsing reaction from section '" + line + "':\n" +
                                         std::string(exc.what()));
            }
        }
        catch (std::exception &exc)
        {
            throw std::runtime_error("While creating collision '" + line + "':\n" + std::string(exc.what()));
        }
    }
}
void EedfGasMixture::loadCollisions(const std::vector<std::string> &files, Grid *energyGrid, bool isExtra)
{
    Log<Message>::Notify("Started loading collisions.");

#ifdef EMSCRIPTEN
    const std::string inputPath{""};
#else
    const std::string inputPath{"../Input/"};
#endif

    for (const std::string &file : files)
    {
        if (file.size() >= 5 && file.substr(file.size() - 5) == ".json")
        {
            const json_type cnf = read_json_from_file(inputPath + file);
            // 1. read the states
            const json_type &scnf = cnf.at("states");
            for (json_type::const_iterator it = scnf.begin(); it != scnf.end(); ++it)
            {
                const StateEntry entry{entryFromJSON(*it)};
                ensureState(entry);
            }
            // 2. read the processes
            for (json_type::const_iterator it = cnf.at("processes").begin(); it != cnf.at("processes").end(); ++it)
            {
                createCollision(*it, energyGrid, isExtra);
            }
        }
        else
        {
            loadCollisions(inputPath + file, energyGrid, isExtra);
        }
    }
#if 0
        std::cout << "List of collisions" << std::endl;
        for (const auto& c : m_collisions)
        {
            std::cout << *c << std::endl;
        }
#endif
    Log<Message>::Notify("Finished loading collisions.");
}

void EedfGasMixture::createCollision(const json_type &pcnf, Grid *energyGrid, bool isExtra)
{
    try
    {
        const json_type &rcnf = pcnf.at("reaction");

        std::vector<uint16_t> lhsCoeffs, rhsCoeffs;

        // 1. Create vectors of pointers to the states that appear
        //    on the left and right-hand sides of the process.
        std::vector<State *> lhsStates;
        std::vector<State *> rhsStates;
        for (const auto &t : rcnf.at("lhs"))
        {
            const std::string &stateName(t.at("state"));
            State *state = findStateById(stateName);
            if (!state)
            {
                throw std::runtime_error("Could not find state '" + stateName + "'.");
            }
            lhsStates.push_back(state);
            lhsCoeffs.push_back(t.contains("count") ? t.at("count").get<int>() : 1);
        }
        for (const auto &t : rcnf.at("rhs"))
        {
            const std::string &stateName(t.at("state"));
            State *state = findStateById(stateName);
            if (!state)
            {
                throw std::runtime_error("Could not find state '" + stateName + "'.");
            }
            rhsStates.push_back(state);
            rhsCoeffs.push_back(t.contains("count") ? t.at("count").get<int>() : 1);
        }

        const CollisionType type = getCollisionTypeFromTypeTagArray(rcnf.at("type_tags"));
        const bool reverseAlso = rcnf.at("reversible");

        // 2. Create the collision object (but do not configure a CrossSection
        //    object yet).
        std::unique_ptr<Collision> collision{
            new Collision(type, lhsStates, lhsCoeffs, rhsStates, rhsCoeffs, reverseAlso)};
        Log<Message>::Notify(*collision);
        // 3. See if we already have a collision of the same type with the same lhs and rhs.
        //    If we do, isue a warning, discard the newly created collision object and return
        //    nullptr. (NOTE: collision is a unique_ptr and will delete the collision object,
        //    since we do not release that before returning from this function..)
        for (const auto &c : m_collisions)
        {
            if (c->is_same_as(*collision))
            {
                Log<DoubleCollision>::Warning(*c);
                collision.release();
                break;
            }
        }
        // 4. Register the collision with the (single) gas that is the target of the
        //    binary electron-heavy process, add it to out own list of collisions
        //    and return the (released) pointer to our caller.
        //    That the left-hand side of the process is indeed of the form 'e + X' is
        //    checked by the EedfCollision constructor, and the target 'X' is returned
        //    by its member getTarget().
        /** \todo Eliminate the following cast at some moment: we can do
         *        EedfGas& eedfGas = *ensureGas(collision->getTarget()->gas().name);
         *        but could check/assert that the pointer is the same.
         */
        if (collision.get())
        {
            EedfGas &eedfGas = static_cast<EedfGas &>(collision->getTarget()->gas());
            eedfGas.addCollision(collision.get(), isExtra);
            m_collisions.push_back(collision.get());
            hasCollisions[static_cast<uint8_t>(collision->type)] = true;

            const bool isElasticOrEffective =
                (collision->type == CollisionType::effective || collision->type == CollisionType::elastic);
            collision->crossSection.reset(new CrossSection(energyGrid, isElasticOrEffective, pcnf));
            collision.release();
        }
    }
    catch (std::exception &exc)
    {
        throw std::runtime_error("Error while parsing reaction from section '" + pcnf.dump(1) + "':\n" +
                                 std::string(exc.what()));
    }
}

void EedfGasMixture::loadGasProperties(const GasPropertiesSetup &setup)
{
    GasMixture::loadGasProperties(setup);
    readGasPropertyFile(gases(), setup.OPBParameter, "OPBParameter", false,
        [](EedfGas& gas, double value) { gas.OPBParameter=value; } );
}

void EedfGasMixture::loadGasProperties(const json_type &cnf)
{
    GasMixture::loadGasProperties(cnf);
    readGasProperty(gases(), cnf, "OPBParameter", false,
        [](EedfGas& gas, double value) { gas.OPBParameter=value; } );
}

void EedfGasMixture::evaluateTotalAndElasticCS()
{
    elasticCrossSection.setZero(grid->nCells() + 1);
    totalCrossSection.setZero(grid->nCells() + 1);

    for (auto &gas : gases())
    {
        if (gas->isDummy())
            continue;

        double massRatio = Constant::electronMass / gas->mass;
        for (auto &collision : gas->collisions[static_cast<uint8_t>(CollisionType::elastic)])
        {
            elasticCrossSection += *collision->crossSection * (collision->getTarget()->density * massRatio);
        }

        for (auto i = static_cast<uint8_t>(CollisionType::elastic); i < gas->collisions.size(); ++i)
        {
            for (auto &collision : gas->collisions[i])
            {
                totalCrossSection += *collision->crossSection * collision->getTarget()->density;

                if (collision->isReverse)
                {
                    Vector superElastic;
                    collision->superElastic(grid->getNodes(), superElastic);

                    totalCrossSection += superElastic * collision->m_rhsHeavyStates[0]->density;
                }
            }
        }
    }
}

void EedfGasMixture::addCARGases(const std::vector<std::string> &CARVector)
{
    for (const auto &gasName : CARVector)
    {
        EedfGas *gas = this->findGas(gasName);

        if (gas == nullptr)
        {
            Log<CARForNonExistent>::Warning(gasName);
            continue;
        }

        gas->checkCARConditions();

        CARGases.emplace_back(gas);
    }
}

/** \bug This ADDS rate coefficients to the previously calulated ones.
 */
void EedfGasMixture::evaluateRateCoefficients(const Vector &eedf)
{
    for (auto &gas : gases())
    {
        for (auto &collVec : gas->collisions)
        {
            for (auto &collision : collVec)
            {
                if (collision->crossSection->threshold() > grid->uMax())
                    continue;
                rateCoefficients.emplace_back(collision->evaluateRateCoefficient(eedf));
            }
        }
        for (auto &collVec : gas->collisionsExtra)
        {
            for (auto &collision : collVec)
            {

                if (collision->crossSection->threshold() > grid->uMax())
                    continue;
                rateCoefficientsExtra.emplace_back(collision->evaluateRateCoefficient(eedf));
            }
        }
    }
}

} // namespace loki
