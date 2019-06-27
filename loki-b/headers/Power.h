//
// Created by daan on 24-06-2019.
//

#ifndef LOKI_CPP_POWER_H
#define LOKI_CPP_POWER_H

#include "Enumeration.h"

struct CollPower {
    double ine{0.}, sup{0.};

    void operator+=(const CollPower &other) {
        ine += other.ine;
        sup += other.sup;
    }
};

struct GasPower {
    double excitationIne{0.}, excitationSup{0.},
            excitationNet{0.}, vibrationalIne{0.}, vibrationalSup{0.},
            vibrationalNet{0.}, rotationalIne{0.}, rotationalSup{0.},
            rotationalNet{0.}, ionizationIne{0.}, attachmentIne{0.};
};

struct Power {
    double field{0.}, elasticNet{0.}, elasticGain{0.},
            elasticLoss{0.}, carNet{0.}, carGain{0.},
            carLoss{0.}, excitationIne{0.}, excitationSup{0.},
            excitationNet{0.}, vibrationalIne{0.}, vibrationalSup{0.},
            vibrationalNet{0.}, rotationalIne{0.}, rotationalSup{0.},
            rotationalNet{0.}, ionizationIne{0.}, attachmentIne{0.},
            inelastic{0.}, superelastic{0.}, eDensGrowth{0.},
            electronElectron{0.}, balance, relativeBalance, reference;

    void operator+=(const GasPower &gasPower) {
        // Apply fancy trick to iterate over member variables of a struct (given they are of the same type).
        auto *itThis = (double *)this;
        auto *it = (double *)&gasPower;

        for (uint32_t i = 0; i < 11; ++i) {
            itThis[i+7] += it[i];
        }
    }
};


#endif //LOKI_CPP_POWER_H
