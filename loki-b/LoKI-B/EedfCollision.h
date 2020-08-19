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

    //class EedfCollision : public Collision<State<EedfGas>>
    class EedfCollision : public Collision<GasBase::State>
    {

        // The raw cross section data and threshold is stored in
        // the CrossSection object

    public:
        using EedfState = GasBase::State;

        // DONE: Inelastic and superelastic rate coefficient variables should be here
        double ineRateCoeff{0.};
        double supRateCoeff{0.};
        std::unique_ptr<CrossSection> crossSection;

        // TODO: Find out the most effective way to pass vectors to this constructor and then to the base class.
        //  Should we use r-value references and move semantics?
        EedfCollision(Enumeration::CollisionType type, std::vector<EedfState *> &reactants,
                      std::vector<EedfState *> &products, std::vector<uint16_t> &stoiCoeff, bool isReverse);

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
