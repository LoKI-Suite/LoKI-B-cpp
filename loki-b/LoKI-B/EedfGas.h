//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_EEDFGAS_H
#define LOKI_CPP_EEDFGAS_H

#include "LoKI-B/CrossSection.h"
#include "LoKI-B/Enumeration.h"
#include "LoKI-B/GasBase.h"
#include "LoKI-B/Grid.h"
#include "LoKI-B/Power.h"

#include <map>
#include <vector>

// TODO: Allow loading of effective populations from a file.
// TODO: comment EedfGas class

namespace loki
{

class EedfCollision;

class EedfGas : public GasBase
{
  public:
    using EedfState = StateBase;
    using CollisionVector = std::vector<std::unique_ptr<EedfCollision>>;
    // We need to store the collisions per Gas since we need to calculate
    // the mass ratio when evaluating the total and elastic cross-sections.
    std::vector<CollisionVector> collisions, collisionsExtra;
    std::map<EedfState *, double> effectivePopulations;
    double OPBParameter = 0.;

    std::map<EedfState *, std::vector<EedfCollision *>> m_state_collisions;
    std::map<EedfState *, std::vector<EedfCollision *>> m_state_collisionsExtra;

    explicit EedfGas(const std::string &name);

    ~EedfGas();

    void addCollision(EedfCollision *collision, bool isExtra);

    void checkElasticCollisions(State *electron, Grid *energyGrid);

    void checkCARConditions();

    bool isDummy();

    const GasPower &getPower() const;

    void evaluatePower(const IonizationOperatorType ionType, const Vector &eedf) const;

  private:
    GasPower power;

    std::vector<EedfState *> findStatesToUpdate();

    CrossSection *elasticCrossSectionFromEffective(Grid *energyGrid);

    void setDefaultEffPop(EedfState *ground);

    CollPower evaluateConservativePower(const CollisionVector &collisionVector, const Vector &eedf) const;

    CollPower evaluateNonConservativePower(const CollisionVector &collisionVector, const IonizationOperatorType ionType,
                                           const Vector &eedf) const;
};
} // namespace loki

#endif // LOKI_CPP_EEDFGAS_H
