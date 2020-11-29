//
// Created by daan on 21-5-19.
//

#ifndef LOKI_CPP_EEDFCOLLISION_H
#define LOKI_CPP_EEDFCOLLISION_H

#include "LoKI-B/Collision.h"
#include "LoKI-B/CrossSection.h"
#include "LoKI-B/EedfGas.h"
#include "LoKI-B/Enumeration.h"
#include "LoKI-B/MacroscopicQuantities.h"
#include "LoKI-B/Power.h"
#include <memory>

namespace loki
{

class EedfCollision : public Collision
{
public:
    using EedfState = GasBase::State;

    EedfCollision(CollisionType type, const StateVector &lhsStates, const CoeffVector &lhsCoeffs,
                  const StateVector &rhsStates, const CoeffVector &rhsCoeffs, bool isReverse);
    ~EedfCollision();
    const EedfState *getTarget() const;
    EedfState *getTarget();
    void superElastic(const Vector &energyData, Vector &result) const;
    PowerTerm evaluateConservativePower(const Vector &eedf) const;
    PowerTerm evaluateNonConservativePower(const Vector &eedf, const IonizationOperatorType ionizationOperatorType,
                                           const double OPBParameter) const;
    /// \todo non-constant because ineRateCoeff, supRateCoeff are modified
    RateCoefficient evaluateRateCoefficient(const Vector &eedf);
    std::string typeAsString() const;
    friend std::ostream &operator<<(std::ostream &os, const EedfCollision &collision);

    double ineRateCoeff() const { return m_ineRateCoeff; }
    double supRateCoeff() const { return m_supRateCoeff; }
private:
    // The raw cross section data and threshold is stored in
    // the CrossSection object
    const StateVector m_lhsHeavyStates;
    const CoeffVector m_lhsHeavyCoeffs;
public:
    /// \todo Make private
    StateVector m_rhsHeavyStates;
private:
    CoeffVector m_rhsHeavyCoeffs;
    double m_ineRateCoeff;
    double m_supRateCoeff;
public:
    std::unique_ptr<CrossSection> crossSection;

};

} // namespace loki

#endif // LOKI_CPP_EEDFCOLLISION_H
