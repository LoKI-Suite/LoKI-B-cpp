//
// Created by daan on 21-5-19.
//

#ifndef LOKI_CPP_EEDFCOLLISION_H
#define LOKI_CPP_EEDFCOLLISION_H

#include "LoKI-B/Collision.h"
#include "LoKI-B/Enumeration.h"
#include "LoKI-B/CrossSection.h"
#include "LoKI-B/Power.h"
#include "LoKI-B/MacroscopicQuantities.h"
#include "LoKI-B/EedfGas.h"
#include <memory>

namespace loki {

    class EedfCollision : public Collision
    {

        // The raw cross section data and threshold is stored in
        // the CrossSection object

        StateVector m_lhsHeavyStates;
    public:
        CoeffVector m_lhsHeavyCoeffs;
        StateVector m_rhsHeavyStates;
        CoeffVector m_rhsHeavyCoeffs;

        using EedfState = GasBase::State;

        // DONE: Inelastic and superelastic rate coefficient variables should be here
        double ineRateCoeff{0.};
        double supRateCoeff{0.};
        std::unique_ptr<CrossSection> crossSection;

        EedfCollision(CollisionType type,
                      const StateVector& lhsStates,
                      const CoeffVector& lhsCoeffs,
                      const StateVector& rhsStates,
                      const CoeffVector& rhsCoeffs,
                      bool isReverse);

        ~EedfCollision();

        const EedfState *getTarget() const;
        EedfState *getTarget();

        void superElastic(const Vector &energyData, Vector &result) const;
        bool is_same_as(CollisionType other_type, const State* other_target, const StateVector& other_products) const;

        friend std::ostream &operator<<(std::ostream &os, const EedfCollision &collision);

        CollPower evaluateConservativePower(const Vector &eedf);

        CollPower evaluateNonConservativePower(const Vector &eedf, const IonizationOperatorType ionizationOperatorType,
                                               const double OPBParameter);

        RateCoefficient evaluateRateCoefficient(const Vector &eedf);

        std::string typeAsString() const;
    };
}


#endif //LOKI_CPP_EEDFCOLLISION_H
