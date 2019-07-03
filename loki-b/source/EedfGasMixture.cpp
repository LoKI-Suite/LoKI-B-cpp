//
// Created by daan on 20-5-19.
//

#include "EedfGasMixture.h"

// DONE: productStoiCoeff is only for the products
// DONE: when parsing an LXCat Section first look for the PARAM keyword (since elastic collisions will not have
//  a threshold).

namespace loki {
    using namespace Enumeration;

    EedfGasMixture::EedfGasMixture(Grid *grid) : grid(grid) {}

    void EedfGasMixture::initialize(const ElectronKineticsSetup &setup, const WorkingConditions *workingConditions) {

        loadCollisions(setup.LXCatFiles, grid);

        if (!setup.LXCatFilesExtra.empty())
            loadCollisions(setup.LXCatFilesExtra, grid, true);

        // DONE: Load Gas properties
        this->loadGasProperties(setup.gasProperties);

        // DONE: Load State properties
        this->loadStateProperties(setup.stateProperties, workingConditions);

        this->evaluateStateDensities();

        for (auto *gas : gasses) {
            if (!gas->isDummy())
                gas->checkElasticCollisions(grid);
        }

        this->addCARGasses(setup.CARgases);

//        this->evaluateTotalAndElasticCS();
    }

    void EedfGasMixture::loadCollisions(const std::vector<std::string> &files, Grid *energyGrid, bool isExtra) {
        const std::string inputPath{"../Input/"};

        const std::regex reParam(R"(PARAM\.:)");
        const std::regex reThreshold(R"(E = (\d*\.?\d*) eV)");
        const std::regex reProcess(R"(\[(.+?)(<->|->)(.+?), (\w+)\])");

        for (const std::string &file : files) {
            std::ifstream in(inputPath + file);

            if (!in.is_open()) {
                Log<FileError>::Warning(file);
                continue;
            }

            std::string line, threshold;
            std::smatch mThreshold, mProcess;

            while (std::getline(in, line)) {
                if (!std::regex_search(line, reParam))
                    continue;

                if (std::regex_search(line, mThreshold, reThreshold)) {
                    threshold = mThreshold[1];
                } else {
                    threshold = "0.";
                }

                if (!std::getline(in, line) ||
                    !std::regex_search(line, mProcess, reProcess)) {
                    Log<LXCatError>::Error(file);
                }

                CollisionEntry collisionEntry = parseLXCatEntry({mProcess[1], mProcess[2],
                                                                 mProcess[3], mProcess[4],
                                                                 threshold});

                auto *collision = createCollision(collisionEntry);

                if (linkCollision(collision, isExtra)) {
                    collision->crossSection = new CrossSection(collisionEntry.threshold, energyGrid, in);
                }
            }
        }
    }

    CollisionEntry EedfGasMixture::parseLXCatEntry(RawLXCatEntry &&entry) {
        CollisionEntry collisionEntry{};

        if (!Parse::entriesFromString(entry.reactants, collisionEntry.reactants))
            Log<LXCatError>::Error(entry.reactants);
        if (!Parse::entriesFromString(entry.products, collisionEntry.products, &collisionEntry.stoiCoeff))
            Log<LXCatError>::Error(entry.products);

        collisionEntry.type = Parse::collisionTypeFromString(entry.type);

        std::stringstream ss(entry.threshold);
        if (!(ss >> collisionEntry.threshold))
            Log<Message>::Error("Non numerical threshold value.");

        collisionEntry.isReverse = (entry.isReverse[0] == '<');

        return collisionEntry;
    }

    bool EedfGasMixture::linkCollision(EedfCollision *collision, bool isExtra) {
        // Linking the newly created collision to the relevant states and gasses

        auto *target = collision->getTarget();

        Log<Message>::Notify(*collision);

        // TODO: threshold <= umax check should probably be performed here.

        // Check if the collision already exists. If so delete the new entry and return false.
        if (target->hasCollision(collision, isExtra)) {
            Log<DoubleCollision>::Warning(*collision);
            delete collision;
            return false;
        }

        target->addCollision(collision, isExtra);
        target->gas->addCollision(collision, isExtra);

        return true;
    }

    void EedfGasMixture::loadGasProperties(const GasPropertiesSetup &setup) {
        std::string fileBuffer;
        GAS_PROPERTY(OPBParameter)

        GasMixture::loadGasProperties(setup);
    }

    void EedfGasMixture::evaluateTotalAndElasticCS() {
        elasticCrossSection.setZero(grid->cellNumber + 1);
        totalCrossSection.setZero(grid->cellNumber + 1);

        for (auto *gas : gasses) {
            if (gas->isDummy()) continue;

            double massRatio = Constant::electronMass / gas->mass;

            for (auto *collision : gas->collisions[(uint8_t) CollisionType::elastic]) {
                elasticCrossSection += *collision->crossSection * (collision->getTarget()->density * massRatio);
            }

            for (auto i = (uint8_t) CollisionType::elastic; i < gas->collisions.size(); ++i) {
                for (auto *collision : gas->collisions[i]) {
                    totalCrossSection += *collision->crossSection * collision->getTarget()->density;

                    if (collision->isReverse) {
                        Vector superElastic;
                        collision->superElastic(grid->getNodes(), superElastic);

                        totalCrossSection += superElastic * collision->products[0]->density;
                    }
                }
            }
        }
    }

    void EedfGasMixture::addCARGasses(const std::vector<std::string> &CARVector) {
        for (const auto &gasName : CARVector) {
            EedfGas * gas = this->findGas(gasName);

            if (gas == nullptr) {
                Log<CARForNonExistent>::Warning(gasName);
                continue;
            }

            gas->checkCARConditions();

            CARGasses.emplace_back(gas);
        }
    }

    void EedfGasMixture::evaluateRateCoefficients(const Vector &eedf) {
        for (auto *gas : gasses) {
            for (auto &collVec : gas->collisions) {
                for (auto *collision : collVec) {

                    if (collision->crossSection->threshold > grid->getNode(grid->cellNumber))
                        continue;

                    rateCoefficients.emplace_back(collision->evaluateRateCoefficient(eedf));
                }
            }

            for (auto &collVec : gas->extraCollisions) {
                for (auto *collision : collVec) {

                    if (collision->crossSection->threshold > grid->getNode(grid->cellNumber))
                        continue;

                    rateCoefficientsExtra.emplace_back(collision->evaluateRateCoefficient(eedf));
                }
            }
        }
    }
}