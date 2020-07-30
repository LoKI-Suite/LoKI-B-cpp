//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_EEDFGAS_H
#define LOKI_CPP_EEDFGAS_H

#include "Enumeration.h"
#include "Gas.h"
#include "EedfCollision.h"
#include "CrossSection.h"
#include "Grid.h"
#include "Power.h"

#include <vector>
#include <map>

// TODO: Allow loading of effective populations from a file.
// TODO: comment EedfGas class

namespace loki {
    class EedfGas : public Gas<Boltzmann> {
    public:
        using CollisionVector = std::vector<std::unique_ptr<EedfCollision>>;
        // We need to store the collisions per Gas since we need to calculate
        // the mass ratio when evaluating the total and elastic cross-sections.
        std::vector<CollisionVector> collisions, extraCollisions;
        std::map<EedfState *, double> effectivePopulations;
        double OPBParameter = 0.;

        explicit EedfGas(const std::string &name);

        ~EedfGas();

        void addCollision(EedfCollision *collision, bool isExtra);

        void checkElasticCollisions(Grid *energyGrid);

        void checkCARConditions();

        bool isDummy();

        const GasPower &getPower() const;

        void evaluatePower(const IonizationOperatorType ionType, const Vector &eedf) const;

    private:
        GasPower power;

        void findStatesToUpdate(const std::vector<EedfState *> &stateStructure,
                                std::vector<EedfState *> &statesToUpdate);

        CrossSection *elasticCrossSectionFromEffective(Grid *energyGrid);

        void setDefaultEffPop(EedfState *ground);

        CollPower evaluateConservativePower(const CollisionVector& collisionVector, const Vector &eedf) const;

        CollPower evaluateNonConservativePower(const CollisionVector& collisionVector,
                                               const IonizationOperatorType ionType, const Vector &eedf) const;
    };
}


#endif //LOKI_CPP_EEDFGAS_H
