//
// Created by daan on 20-5-19.
//

#ifndef LOKI_CPP_EEDFGASMIXTURE_H
#define LOKI_CPP_EEDFGASMIXTURE_H

#include "LoKI-B/EedfCollision.h"
#include "LoKI-B/EedfState.h"
#include "LoKI-B/EedfGas.h"
#include "LoKI-B/GasMixture.h"
#include "LoKI-B/WorkingConditions.h"
#include "LoKI-B/MacroscopicQuantities.h"
#include "LoKI-B/Traits.h"
#include "LoKI-B/json.h"

#include <vector>

namespace loki {

    class EedfGasMixture : public GasMixture<Boltzmann>
    {
    public:

        explicit EedfGasMixture(Grid *grid);
        /** Initializes the gas mixture by loading the desired collisions from LXCat files.
         *  These files are read from the electron kinetics setup structure. It also
         *  requires a pointer to the energy grid in order to properly initialize the
         *  cross sections of the collisions.
         */
        void initialize(const ElectronKineticsSetup &setup, const WorkingConditions *workingConditions);
        void initialize(const json_type &cnf, const WorkingConditions *workingConditions);

        Vector elasticCrossSection, totalCrossSection;

        Grid *grid{nullptr};

        std::vector<EedfGas *> CARGases;

        std::vector<RateCoefficient> rateCoefficients, rateCoefficientsExtra;

        bool hasCollisions[static_cast<uint8_t>(CollisionType::size)]{false};


        // TODO: comment evaluateTotalAndElasticCS

        void evaluateTotalAndElasticCS();

        void evaluateRateCoefficients(const Vector &eedf);

    private:

        /* -- loadCollisions --
         * Loads the collisions from the LXCat file that is provided as first argument.
         * Furthermore, it needs a pointer to the energy grid and a boolean to indicate
         * whether the collisions are extra, for correct initialization and storage of
         * the collisions.
         */

        void loadCollisions(const std::string& file, Grid *energyGrid, bool isExtra);

        /* -- loadCollisionsJSON --
         * Loads the collisions from the JSON file that is provided as first argument.
         * Furthermore, it needs a pointer to the energy grid and a boolean to indicate
         * whether the collisions are extra, for correct initialization and storage of
         * the collisions.
         */
        void loadCollisionsJSON(const json_type& file, Grid *energyGrid, bool isExtra);

        /* -- loadCollisions --
         * Loads the collisions from files, supplied through a vector of strings that hold
         * the filenames. Furthermore, it needs a pointer to the energy grid and a boolean to
         * indicate whether the collisions are extra, for correct initialization and storage of
         * the collisions. When the file extension is ".json", a JSON object is created from
         * the file and the handling of this file is delegated to member loadCollisionsJSON,
         * for other file types, the legacy LXCat file format is assumed and member
         * loadCollisionsText is called.
         */

        void loadCollisions(const std::vector<std::string> &files, Grid *energyGrid, bool isExtra = false);

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
        void loadGasProperties(const json_type &cnf) override;

        // TODO: comment addCARGases

        void addCARGases(const std::vector<std::string> &CARVector);
    };
}


#endif //LOKI_CPP_EEDFGASMIXTURE_H
