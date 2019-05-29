//
// Created by daan on 20-5-19.
//

#include "EedfGasMixture.h"

// DONE: productStoiCoeff is only for the products
// DONE: when parsing an LXCat Section first look for the PARAM keyword (since elastic collisions will not have
//  a threshold).

namespace loki {
    void EedfGasMixture::initialize(const ElectronKineticsSetup &setup, Grid *energyGrid,
            const WorkingConditions *workingConditions) {

        loadCollisions(setup.LXCatFiles, energyGrid);

        if (!setup.LXCatFilesExtra.empty())
            loadCollisions(setup.LXCatFilesExtra, energyGrid, true);

        // DONE: Load Gas properties
        this->loadGasProperties(setup.gasProperties);

        // TODO: Load State properties
        this->loadStateProperties(setup.stateProperties, workingConditions);
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
}