//
// Created by daan on 20-5-19.
//

#include <stdexcept>
#include <iostream>
#include "LoKI-B/EedfGasMixture.h"
#include "LoKI-B/Constant.h"
#include "LoKI-B/json.h"

// DONE: productStoiCoeff is only for the products
// DONE: when parsing an LXCat Section first look for the PARAM keyword (since elastic collisions will not have
//  a threshold).

namespace loki {

    using namespace Enumeration;

    EedfGasMixture::EedfGasMixture(Grid *grid)
    : grid(grid)
    {
    }
    void EedfGasMixture::initialize(const ElectronKineticsSetup &setup, const WorkingConditions *workingConditions)
    {

        loadCollisions(setup.LXCatFiles, grid);

        if (!setup.LXCatFilesExtra.empty())
            loadCollisions(setup.LXCatFilesExtra, grid, true);

        // DONE: Load Gas properties
        this->loadGasProperties(setup.gasProperties);

        // DONE: Load State properties
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
    void EedfGasMixture::initialize(const json_type &cnf, const WorkingConditions *workingConditions)
    {

        loadCollisions(cnf.at("LXCatFiles"), grid);

        if (cnf.contains("LXCatFilesExtra"))
            loadCollisions(cnf.at("LXCatFilesExtra").get<std::vector<std::string>>(), grid, true);

        // DONE: Load Gas properties
        this->loadGasProperties(cnf.at("gasProperties"));

        // DONE: Load State properties
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
            /** \todo try to make this more exception-safe. The 'new' that is part of
             *  createCollision is quite detached from the delete in linkCollision
             *  (which should be called if something goes wrong on the way).
             */
            auto *collision = createCollision(rcnf);
            const bool isElasticOrEffective = (collision->type == CollisionType::effective ||
                                               collision->type == CollisionType::elastic);
            if (linkCollision(collision, isExtra))
            {
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
            auto *collision = createCollision(mProcess[1], mProcess[2],mProcess[3], mProcess[4]);

            const bool isElasticOrEffective = (collision->type == CollisionType::effective ||
                                               collision->type == CollisionType::elastic);
            if (linkCollision(collision, isExtra))
            {
                collision->crossSection.reset(new CrossSection(threshold, energyGrid,
                                                           isElasticOrEffective, in));
                hasCollisions[static_cast<uint8_t>(collision->type)] = true;
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

    bool EedfGasMixture::linkCollision(EedfCollision *collision, bool isExtra)
    {
        // Linking the newly created collision to the relevant states and gases

        auto *target = collision->getTarget();

        Log<Message>::Notify(*collision);

        // Check if the collision already exists. If so delete the new entry and return false.
        if (target->hasCollision(collision, isExtra))
        {
            Log<DoubleCollision>::Warning(*collision);
            delete collision;
            return false;
        }

        target->addCollision(collision, isExtra);
        target->gas().addCollision(collision, isExtra);

        return true;
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

                        totalCrossSection += superElastic * collision->products[0]->density;
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
            for (auto &collVec : gas->extraCollisions)
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
}
