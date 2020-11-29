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

    explicit EedfGas(const std::string &name);
    ~EedfGas();

    void addCollision(EedfCollision *collision, bool isExtra);
    void checkElasticCollisions(State *electron, Grid *energyGrid);
    bool isDummy() const;
    const GasPower &getPower() const;
    /** \todo Non-constant because m_power is changed. See if m_power must managed here.
     *  The name updatePower would make more clear that this modifies the object.
     */
    const GasPower& evaluatePower(const IonizationOperatorType ionType, const Vector &eedf);
    double OPBParameter() const { return m_OPBParameter; }
    void setOPBParameter(double value) { m_OPBParameter = value; }

    // We need to store the collisions per Gas since we need to calculate
    // the mass ratio when evaluating the total and elastic cross-sections.
    const CollisionVector& collisions(CollisionType type) const { return m_collisions[static_cast<uint8_t>(type)]; }
    const std::vector<CollisionVector>& collisions() const { return m_collisions; }
    const std::vector<CollisionVector>& collisionsExtra() const { return m_collisionsExtra; }
private:
    std::vector<CollisionVector> m_collisions;
    std::vector<CollisionVector> m_collisionsExtra;
    std::map<const EedfState *, double> m_effectivePopulations;
    double m_OPBParameter;
    std::map<const EedfState *, std::vector<EedfCollision *>> m_state_collisions;
    std::map<const EedfState *, std::vector<EedfCollision *>> m_state_collisionsExtra;
    GasPower m_power;

    std::vector<EedfState *> findStatesToUpdate();
    CrossSection *elasticCrossSectionFromEffective(Grid *energyGrid);
    void setDefaultEffPop(EedfState *ground);
    PowerTerm evaluateConservativePower(const CollisionVector &collisionVector, const Vector &eedf) const;
    PowerTerm evaluateNonConservativePower(const CollisionVector &collisionVector, const IonizationOperatorType ionType,
                                           const Vector &eedf) const;
};
} // namespace loki

#endif // LOKI_CPP_EEDFGAS_H
