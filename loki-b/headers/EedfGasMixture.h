//
// Created by daan on 20-5-19.
//

#ifndef LOKI_CPP_EEDFGASMIXTURE_H
#define LOKI_CPP_EEDFGASMIXTURE_H

#include "GasMixture.h"
#include "WorkingConditions.h"
#include "Traits.h"

namespace loki {
    class EedfGasMixture : public GasMixture<Boltzmann> {
    public:

        /* -- initialize --
         * Initializes the gas mixture by loading the desired collisions from LXCat files.
         * These files are read from the electron kinetics setup structure. It also
         * requires a pointer to the energy grid in order to properly initialize the
         * cross sections of the collisions.
         */

        void initialize(const ElectronKineticsSetup &setup, Grid *energyGrid,
                const WorkingConditions *workingConditions);

    private:

        /* -- loadCollisions --
         * Loads the collisions from LXCat files, supplied through a vector of strings that hold
         * the filenames. Furthermore, it needs a pointer to the energy grid and a boolean to
         * indicate whether the collisions are extra, for correct initialization and storage of
         * the collisions.
         */

        void loadCollisions(const std::vector<std::string> &files, Grid *energyGrid, bool isExtra = false);

        /* -- parseLXCatEntry --
         * Parses a RawLXCatEntry into StateEntries which build up a CollisionEntry object. This
         * object acts as a blueprint for a collision.
         */

        CollisionEntry parseLXCatEntry(RawLXCatEntry &&entry);

        /* -- linkCollision --
         * Accepts a pointer to a newly created EedfCollision to the correct vectors of the
         * target gas and states. It also checks whether the collision already exists, if
         * this is the case then the new collision is deleted and the function returns false.
         */

        bool linkCollision(EedfCollision *collision, bool isExtra);

        /* -- loadGasProperties --
         * EedfGas introduces one extra property that has to be set from a file: OPBParameter.
         * This overload sets this parameter and then calls Gas::loadGasProperties to set the
         * rest of its properties.
         */

        void loadGasProperties(const GasPropertiesSetup &setup) override;
    };
}


#endif //LOKI_CPP_EEDFGASMIXTURE_H
