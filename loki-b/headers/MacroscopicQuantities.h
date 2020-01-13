//
// Created by daan on 02-07-2019.
//

#ifndef LOKI_CPP_MACROSCOPICQUANTITIES_H
#define LOKI_CPP_MACROSCOPICQUANTITIES_H

namespace loki {
    struct SwarmParameters {
        double redDiffCoeff{0.}, redMobCoeff{0.}, redTownsendCoeff{0.},
                redAttCoeff{0.}, meanEnergy{0.}, characEnergy{0.},
                Te{0.}, driftVelocity{0.};
    };

    // Forward declaration
    class EedfCollision;

    struct RateCoefficient {
        EedfCollision * collision;
        double inelastic, superelastic;
    };
}

#endif //LOKI_CPP_MACROSCOPICQUANTITIES_H
