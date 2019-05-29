//
// Created by daan on 21-5-19.
//

#ifndef LOKI_CPP_EEDFCOLLISION_H
#define LOKI_CPP_EEDFCOLLISION_H

#include "Collision.h"
#include "Enumeration.h"
#include "CrossSection.h"

namespace loki {
    class EedfCollision : Collision<Boltzmann> {

        // The raw cross section data and threshold is stored in
        // the CrossSection object

        // DONE: Inelastic and superelastic rate coefficient variables should be here
        double ineRateCoeff{0.}, supRateCoef{0.};

    public:
        CrossSection * crossSection{nullptr};

        // TODO: Find out the most effective way to pass vectors to this constructor and then to the base class.
        //  Should we use r-value references and move semantics?
        EedfCollision(Enumeration::CollisionType type, std::vector<EedfState *> &reactants,
                std::vector<EedfState *> &products, std::vector<uint16_t> &stoiCoeff, bool isReverse);

        ~EedfCollision();

        EedfState *getTarget();

        bool operator==(const EedfCollision &other);

        friend std::ostream &operator<<(std::ostream& os, const EedfCollision &collision);
    };
}


#endif //LOKI_CPP_EEDFCOLLISION_H
