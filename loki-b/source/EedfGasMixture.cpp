//
// Created by daan on 20-5-19.
//

#include <stdexcept>
#include <iostream>
#include "LoKI-B/EedfGasMixture.h"
#include "LoKI-B/Constant.h"
#include "LoKI-B/json.h"

namespace loki {

    EedfGasMixture::EedfGasMixture(Grid *grid, const ElectronKineticsSetup &setup, const WorkingConditions *workingConditions)
    : grid(grid)
    {

        loadCollisions(setup.LXCatFiles, grid);

        if (!setup.LXCatFilesExtra.empty())
            loadCollisions(setup.LXCatFilesExtra, grid, true);

        this->loadGasProperties(setup.gasProperties);
        this->loadStateProperties(setup.stateProperties, workingConditions);
        this->evaluateStateDensities();

        for (auto& gas : gases())
        {
            if (!gas->isDummy())
                gas->checkElasticCollisions(grid);
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
        this->loadStateProperties(cnf.at("stateProperties"), workingConditions);
        this->evaluateStateDensities();

        for (auto& gas : gases())
        {
            if (!gas->isDummy())
                gas->checkElasticCollisions(grid);
        }
        if (cnf.contains("CARgases"))
        {
            this->addCARGases(cnf.at("CARgases"));
        }

//        this->evaluateTotalAndElasticCS();
    }

    void EedfGasMixture::loadCollisionsJSON(const json_type& cnf, Grid *energyGrid, bool isExtra)
    {
        unsigned ndx=0;
        for (json_type::const_iterator it = cnf.begin(); it != cnf.end(); ++it)
        {
            const json_type& rcnf = it->at("reaction");
            auto *collision = createCollision(rcnf,isExtra);
            if (collision)
            {
                const bool isElasticOrEffective = (collision->type == CollisionType::effective ||
                                               collision->type == CollisionType::elastic);
                const double threshold = it->contains("threshold") ? it->at("threshold").get<double>() : 0.0;
                collision->crossSection.reset(new CrossSection(threshold, energyGrid,
                                                           isElasticOrEffective, *it));
                hasCollisions[static_cast<uint8_t>(collision->type)] = true;
            }

        }
    }
    void EedfGasMixture::loadCollisions(const std::string& file, Grid *energyGrid, bool isExtra)
    {
        const std::regex reParam(R"(PARAM\.:)");
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
         *  1. Read the arguments of the parameter line, see if a threshold is specified ("E =  <double> eV"),
         *     Set the threshold to 0 otherwise.
         *  2. The next line should describe the process: it must be of the form [<lhs> -> <rhs>, <type>].
         *  3. A Collision is created from the process description information
         *  4. (Maybe) a cross section object is attached to the collision object, based on subsequent lines.
         *
         */
        while (std::getline(in, line))
        {
        try {
            if (!std::regex_search(line, reParam))
                continue;
            // todo: warn if a threshold is specified, but isExtra==false? Then it will not be used, it seems.
            double threshold = 0.0;
            if (std::regex_search(line, mThreshold, reThreshold))
            {
                std::stringstream ss(mThreshold[1]);
                if (!(ss >> threshold))
                    Log<Message>::Error("Non numerical threshold value.");
            }
            if (!std::getline(in, line) ||
                !std::regex_search(line, mProcess, reProcess))
            {
                Log<LXCatError>::Error(file);
            }

            // arguments: smth. like "He + e", "->" or "<->", "He + e", "Elastic"
            auto *collision = createCollision(mProcess[1], mProcess[2],mProcess[3], mProcess[4],isExtra);

            if (collision)
            {
                const bool isElasticOrEffective = (collision->type == CollisionType::effective ||
                                               collision->type == CollisionType::elastic);
                collision->crossSection.reset(new CrossSection(threshold, energyGrid,
                                                           isElasticOrEffective, in));
                hasCollisions[static_cast<uint8_t>(collision->type)] = true;
            }
        }
        catch(std::exception& exc)
        {
            throw std::runtime_error("While creating collision '" + line + "':\n" + std::string(exc.what()));
        }
        }
    }
    void EedfGasMixture::loadCollisions(const std::vector<std::string> &files, Grid *energyGrid, bool isExtra) {

        std::cout << "Starting loading collisions" << std::endl;
        const std::string inputPath{"../Input/"};
        for (const std::string &file : files)
        {
            if (file.size()>=5 && file.substr(file.size()-5)==".json")
            {
                json_type cnf = read_json_from_file(inputPath + file);
                loadCollisionsJSON(cnf, energyGrid, isExtra);
            }
            else
            {
                loadCollisions(inputPath + file, energyGrid, isExtra);
            }
        }
        std::cout << "Finished loading collisions" << std::endl;
    }

    EedfGasMixture::Collision* EedfGasMixture::createCollision(
                CollisionType entry_type,
                const std::vector<StateEntry>& entry_lhsStates,
                const std::vector <uint16_t>& entry_lhsCoeffs,
                const std::vector<StateEntry>& entry_rhsStates,
                const std::vector <uint16_t>& entry_rhsCoeffs,
                bool reverse_also,
                bool isExtra)
    {

        std::vector<State*> lhsStates;
        std::vector<State*> rhsStates;
        std::set<GasBase*> targetGases;

        for (auto &stateEntry : entry_lhsStates)
        {
            auto *state = lhsStates.emplace_back(ensureState(stateEntry));
            targetGases.insert(&state->gas());
        }

        if (targetGases.size() != 1)
            Log<Message>::Error("Multiple target gases in a single collision.");

        for (auto &stateEntry : entry_rhsStates)
        {
            rhsStates.emplace_back(ensureState(stateEntry));
        }
        /// \todo Eliminate this cast at some moment
        // reference to the one-and-only target gas for this collision
        // (that there is only one has been established a few lines up).
        EedfGas& eedfGas = static_cast<EedfGas&>(**targetGases.begin());

        /** \todo Check if we have this process already without create a Collision object.
         *        That also allows us to get rid of the EEDFCollision's operator==
         */

        // Linking the newly created collision to the relevant states and gases

        auto *target = lhsStates.front();

        // small change: before there were separate lists (extra / non-extra). Duplication
        // was not detected if the same collision would be added to both lists.
        for (const auto& c : m_collisions)
        {
            if (c->is_same_as(entry_type, target, rhsStates))
            {
                Log<DoubleCollision>::Warning(*c);
                return nullptr;
            }
        }
        std::unique_ptr<Collision> collision{new Collision(entry_type,
                    lhsStates, entry_lhsCoeffs,
                    rhsStates, entry_rhsCoeffs,
                    reverse_also)};
        Log<Message>::Notify(*collision);
        eedfGas.addCollision(collision.get(), isExtra);
        m_collisions.push_back(collision.get());
        return collision.release();
    }

    EedfGasMixture::Collision* EedfGasMixture::createCollision(const std::string& lhs, const std::string& sep, const std::string& rhs, const std::string& type,bool isExtra)
    {
//#define LOKIB_USE_OLD_STATEENTRY_PARSER
#ifdef LOKIB_USE_OLD_STATEENTRY_PARSER
        std::vector <StateEntry> entry_lhsStates, entry_rhsStates;
        std::vector <uint16_t> entry_lhsCoeffs,entry_rhsCoeffs;

        const CollisionType entry_type = getCollisionType(type);
        const bool entry_isReverse = (sep[0] == '<');

        /// \todo Check that the stoichiometric coefficients are indeed 1 (or unspecified), since these are ignored.
        entriesFromStringOld(lhs, entry_lhsStates, &entry_lhsCoeffs);
        entriesFromStringOld(rhs, entry_rhsStates, &entry_rhsCoeffs);
        assert(entry_lhsStates.size()==entry_lhsCoeffs.size());
        if (entry_lhsStates.size() != 1 || entry_lhsCoeffs[0] != 1)
        {
            Log<Message>::Error("Expected one target in collision's left-hand side '" + lhs + "'.");
        }
#else
        // if we use the new parser, also electron states that participate in
        // reactions are recognized and created. Since, at present, the
        // EedfCollision class expects to receive only the heavies on both
        // sides of the reaction, we remove the electrons from the state
        // and coefficient arrays. For the left-hand side this is done in
        // a special way, since we also need to check that the reaction is
        // of the form 'e + X -> ...': it must be a binary encounter of an
        // electron ans a heavy particle.

        /** \todo It seems to make sense to do this stripping in the EedfCollision
         *        constructor, since that class really decide that it only wants
         *        to see the heavies. The Collision base class could then have
         *        the full species/coefficient lists.
         */

        std::vector <StateEntry> entry_lhsStates, entry_rhsStates;
        std::vector <uint16_t> entry_lhsCoeffs,entry_rhsCoeffs;

        const CollisionType entry_type = getCollisionType(type);
        const bool entry_isReverse = (sep[0] == '<');

        entriesFromStringNew(lhs, entry_lhsStates, &entry_lhsCoeffs);
        entriesFromStringNew(rhs, entry_rhsStates, &entry_rhsCoeffs);

        // for the left hand side we expect 'e + X' or 'X + e'.
        assert(entry_lhsStates.size()==entry_lhsCoeffs.size());
        if (entry_lhsStates.size() != 2
            || entry_lhsCoeffs[0]!=1
            || entry_lhsCoeffs[1]!=1
            || (entry_lhsStates[0].gasName=="e" && entry_lhsStates[1].gasName=="e")
            || (entry_lhsStates[0].gasName!="e" && entry_lhsStates[1].gasName!="e")
        )
        {
            Log<Message>::Error("Expected an electron-heavy process on left-hand side '" + lhs + "'.");
        }
        unsigned electron_ndx = entry_lhsStates[0].gasName=="e" ? 0 : 1;
        // remove the electron entry, the present Eedf collision code only wants to see the target
        entry_lhsStates.erase(entry_lhsStates.begin()+electron_ndx);
        entry_lhsCoeffs.erase(entry_lhsCoeffs.begin()+electron_ndx);
        assert(entry_lhsStates.size()==1);
        assert(entry_lhsCoeffs.size()==1);
        // now remove all electrons from the right-hand side (since EedfCollision wants that).
        // beware of iterator invalidation
        for (unsigned i=0; i!= entry_rhsStates.size(); /**/)
        {
            if (entry_rhsStates[i].gasName=="e")
            {
                entry_rhsStates.erase(entry_rhsStates.begin()+i);
                entry_rhsCoeffs.erase(entry_rhsCoeffs.begin()+i);
            }
            else
            {
                ++i;
            }
        }
#endif // LOKIB_USE_OLD_STATEENTRY_PARSER
        return createCollision(entry_type,entry_lhsStates,entry_lhsCoeffs,entry_rhsStates,entry_rhsCoeffs,entry_isReverse,isExtra);
    }

    EedfGasMixture::Collision* EedfGasMixture::createCollision(const json_type& rcnf,bool isExtra)
    {
        std::vector <StateEntry> entry_lhsStates, entry_rhsStates;
        std::vector <uint16_t> entry_lhsCoeffs, entry_rhsCoeffs;

        if (!entriesFromJSON(rcnf.at("lhs"), entry_lhsStates, &entry_lhsCoeffs))
            Log<LXCatError>::Error(rcnf.at("lhs").dump(2));
        if (!entriesFromJSON(rcnf.at("rhs"), entry_rhsStates, &entry_rhsCoeffs))
            Log<LXCatError>::Error(rcnf.at("rhs").dump(2));
        const CollisionType entry_type = getCollisionType(rcnf.at("type"));
        const bool entry_isReverse = rcnf.at("superelastic");

        return createCollision(entry_type,entry_lhsStates,entry_lhsCoeffs,entry_rhsStates,entry_rhsCoeffs,entry_isReverse,isExtra);
    }

    void EedfGasMixture::loadGasProperties(const GasPropertiesSetup &setup)
    {
        std::string fileBuffer;
        GAS_PROPERTY(gases(),OPBParameter)

        GasMixture::loadGasProperties(setup);
    }

    void EedfGasMixture::loadGasProperties(const json_type &cnf)
    {
        std::string fileBuffer;
        GAS_PROPERTY_JSON(gases(),cnf,OPBParameter)

        GasMixture::loadGasProperties(cnf);
    }

    void EedfGasMixture::evaluateTotalAndElasticCS()
    {
        elasticCrossSection.setZero(grid->cellNumber + 1);
        totalCrossSection.setZero(grid->cellNumber + 1);

        for (auto& gas : gases())
        {
            if (gas->isDummy())
              continue;

            double massRatio = Constant::electronMass / gas->mass;
            for (auto& collision : gas->collisions[static_cast<uint8_t>(CollisionType::elastic)])
            {
                elasticCrossSection += *collision->crossSection * (collision->getTarget()->density * massRatio);
            }

            for (auto i = static_cast<uint8_t>(CollisionType::elastic); i < gas->collisions.size(); ++i)
            {
                for (auto& collision : gas->collisions[i])
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

    void EedfGasMixture::evaluateRateCoefficients(const Vector &eedf)
    {
        for (auto& gas : gases())
        {
            for (auto &collVec : gas->collisions)
            {
                for (auto& collision : collVec)
                {
                    if (collision->crossSection->threshold > grid->getNode(grid->cellNumber))
                        continue;
                    rateCoefficients.emplace_back(collision->evaluateRateCoefficient(eedf));
                }
            }
            for (auto &collVec : gas->collisionsExtra)
            {
                for (auto& collision : collVec)
                {

                    if (collision->crossSection->threshold > grid->getNode(grid->cellNumber))
                        continue;
                    rateCoefficientsExtra.emplace_back(collision->evaluateRateCoefficient(eedf));
                }
            }
        }
    }

} // namespace loki
