//
// Created by daan on 20-5-19.
//

#ifndef LOKI_CPP_EEDFGASMIXTURE_H
#define LOKI_CPP_EEDFGASMIXTURE_H

#include "LoKI-B/EedfCollision.h"
#include "LoKI-B/EedfGas.h"
#include "LoKI-B/GasMixture.h"
#include "LoKI-B/WorkingConditions.h"
#include "LoKI-B/MacroscopicQuantities.h"
#include "LoKI-B/json.h"

#include <vector>

namespace loki {

    class EedfGasMixture : public GasMixture<EedfGas>
    {
    public:
        using Collision = EedfCollision;

        /** Initializes the gas mixture by loading the desired collisions from LXCat files.
         *  These files are read from the electron kinetics setup structure. It also
         *  requires a pointer to the energy grid in order to properly initialize the
         *  cross sections of the collisions.
         */
        EedfGasMixture(Grid *grid, const ElectronKineticsSetup &setup, const WorkingConditions *workingConditions);
        EedfGasMixture(Grid *grid, const json_type &cnf, const WorkingConditions *workingConditions);

        Vector elasticCrossSection, totalCrossSection;

        Grid *grid{nullptr};

        std::vector<EedfGas *> CARGases;

        std::vector<RateCoefficient> rateCoefficients, rateCoefficientsExtra;

        bool hasCollisions[static_cast<uint8_t>(CollisionType::size)]{false};


        // TODO: comment evaluateTotalAndElasticCS

        void evaluateTotalAndElasticCS();

        void evaluateRateCoefficients(const Vector &eedf);

    private:

        /** Creates a collision based on a provided CollisionEntry object. The
         *  gases and states involved in the collision are first created and
         *  added to the mixture. Then a Collision object is created and its
         *  pointer is returned.
         */
        // arguments: smth. like "He + e", "->", "He + e", "Elastic"
        Collision* createCollision(const std::string& lhs, const std::string& sep, const std::string& rhs, const std::string& type, bool isExtra);
        Collision* createCollision(const json_type& rcnf, bool isExtra);
        /** Create a collision if it does not already exists, and attaches it to the target gas/state.
         */
        Collision* createCollision(
                Enumeration::CollisionType entry_type,
                std::vector<StateEntry> entry_reactants,
                std::vector<StateEntry> entry_products,
                std::vector <uint16_t> entry_stoiCoeff,
                bool reverse_also,
                bool isExtra);
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

        /* -- loadGasProperties --
         * EedfGas introduces one extra property that has to be set from a file: OPBParameter.
         * This overload sets this parameter and then calls Gas::loadGasProperties to set the
         * rest of its properties.
         */

        void loadGasProperties(const GasPropertiesSetup &setup) override;
        void loadGasProperties(const json_type &cnf) override;

        // TODO: comment addCARGases

        void addCARGases(const std::vector<std::string> &CARVector);

        std::vector<const Collision*> m_collisions;
    };
}


#endif //LOKI_CPP_EEDFGASMIXTURE_H
