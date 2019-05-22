//
// Created by daan on 20-5-19.
//

#ifndef LOKI_CPP_EEDFGASMIXTURE_H
#define LOKI_CPP_EEDFGASMIXTURE_H

#include "GasMixture.h"
#include "Traits.h"

namespace loki {
    class EedfGasMixture : public GasMixture<Boltzmann> {
    public:
        // TODO: Comment initialize
        void initialize(const ElectronKineticsSetup &setup, const Grid *energyGrid);

    private:
        // TODO: Comment loadCollisions
        void loadCollisions(const std::vector<std::string> &files, const Grid *energyGrid, bool isExtra = false);

        /* -- parseLXCatEntry --
         * Creates StateEntries and adds them to a CollisionEntry. The Collision is then created
         * through the CreateCollision function this function will then return a pointer to said
         * collision.
         */

        CollisionEntry parseLXCatEntry(RawLXCatEntry entry, std::ifstream &in);

        /* -- linkCollision --
         * Adds a pointer to the newly created EedfCollision to the correct vectors of the target
         * gas and states. It also checks whether the collision already exists, if this is the
         * case then the new collision is deleted.
         */

        bool linkCollision(EedfCollision *collision, bool isExtra);


    };
}


#endif //LOKI_CPP_EEDFGASMIXTURE_H
