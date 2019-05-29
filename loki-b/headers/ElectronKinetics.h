//
// Created by daan on 13-5-19.
//

#ifndef LOKI_CPP_ELECTRONKINETICS_H
#define LOKI_CPP_ELECTRONKINETICS_H

#include "Setup.h"
#include "EedfGasMixture.h"
#include "Grid.h"
#include "WorkingConditions.h"

namespace loki {
    using namespace Enumeration;

    class ElectronKinetics {
        EedfType eedfType;
        uint8_t shapeParameter;
        IonizationOperatorType ionizationOperatorType;
        GrowthModelType growthModelType;
        bool includeEECollisions;

        const WorkingConditions *workingConditions;

        Grid grid;

        EedfGasMixture mixture;

    public:
        explicit ElectronKinetics(const ElectronKineticsSetup &setup, const WorkingConditions *workingConditions);
        ~ElectronKinetics() = default;

        // Copying this object is not allowed.
        ElectronKinetics(const ElectronKinetics &other) = delete;
    };
}


#endif //LOKI_CPP_ELECTRONKINETICS_H
