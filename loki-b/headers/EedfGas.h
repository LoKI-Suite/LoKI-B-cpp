//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_EEDFGAS_H
#define LOKI_CPP_EEDFGAS_H

#include "Gas.h"
#include "EedfState.h"
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

        // We need to store the collisions per Gas since we need to calculate
        // the mass ratio when evaluating the total and elastic cross-sections.
        std::vector<std::vector<EedfCollision *>> collisions, extraCollisions;
        std::map<EedfState *, double> effectivePopulations;
        double OPBParameter = 0.;

        explicit EedfGas(const std::string &name);

        ~EedfGas();

        void addCollision(EedfCollision *collision, bool isExtra);

        void checkElasticCollisions(Grid *energyGrid);

        void checkCARConditions();

        bool isDummy();

        const GasPower &getPower() const;

        void evaluatePower(const IonizationOperatorType ionType, const Vector &eedf);

    private:
        GasPower power;

        void findStatesToUpdate(const std::vector<EedfState *> &stateStructure,
                                std::vector<EedfState *> &statesToUpdate);

        CrossSection *elasticCrossSectionFromEffective(Grid *energyGrid);

        void setDefaultEffPop(EedfState *ground);

        CollPower evaluateConservativePower(std::vector<EedfCollision *> &collisionVector, const Vector &eedf);

        CollPower evaluateNonConservativePower(std::vector<EedfCollision *> &collisionVector,
                                               const IonizationOperatorType ionType, const Vector &eedf);
    };
}


#endif //LOKI_CPP_EEDFGAS_H
