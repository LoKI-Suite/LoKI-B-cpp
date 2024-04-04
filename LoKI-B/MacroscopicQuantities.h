//
// Created by daan on 02-07-2019.
//

#ifndef LOKI_CPP_MACROSCOPICQUANTITIES_H
#define LOKI_CPP_MACROSCOPICQUANTITIES_H

namespace loki {

struct SwarmParameters
{
    double redDiffCoeff{0.};
    double redMobCoeff{0.};
    double redTownsendCoeff{0.};
    double redAttCoeff{0.};
    double redDiffCoeffEnergy{-1.};
    double redMobilityEnergy{-1.};
    double meanEnergy{0.};
    double characEnergy{0.};
    double Te{0.};
    double driftVelocity{0.};
};

// Forward declaration
class EedfCollision;

struct RateCoefficient
{
    const EedfCollision *collision;
    double inelastic;
    double superelastic;
};

} // namespace loki

#endif // LOKI_CPP_MACROSCOPICQUANTITIES_H
