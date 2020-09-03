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
        /** \todo Update docs, this now does ALL the work...
         *  This member contains the bits of the collision creation code
         *  that are shared by the legacy and JSON set up code. It adds
         *  the states that are mentioned on the LHS and RHS of the equation,
         *  then creates the collision object. If an eqiuvalent object
         *  already exists (same particles on both sides, same type), the
         *  object is discarded and a nullptr is returned. Otherwise,
         *  addCollision is called on the target gas of this collision, the
         *  collision is added to our own list of all collisions, the
         *  hasCollisions entry for the particular type of collision is
         *  set to true, and the collision pointer is returned. In this
         *  case, the caller will still need to configure a cross section
         *  object for this collision. That task is not part of this
         *  function since it depends on the input style legacy/JSON.
         */
        void createCollision(const json_type& pcnf, Grid *energyGrid, bool isExtra);

        /* -- loadCollisions --
         * Loads the collisions from the file that is provided as first argument.
         * Furthermore, it needs a pointer to the energy grid and a boolean to indicate
         * whether the collisions are extra, for correct initialization and storage of
         * the collisions.
         */

        void loadCollisions(const std::string& file, Grid *energyGrid, bool isExtra);

        /* -- loadCollisions --
         * Loads the collisions from files, supplied through a vector of strings that hold
         * the filenames. Furthermore, it needs a pointer to the energy grid and a boolean to
         * indicate whether the collisions are extra, for correct initialization and storage of
         * the collisions. When the file extension is ".json", a JSON object is created from
         * the file and the handling of this file is delegated to member loadCollisionsJSON,
         * for other file types, the legacy LXCat file format is assumed and member
         * loadCollisions is called for the file.
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
