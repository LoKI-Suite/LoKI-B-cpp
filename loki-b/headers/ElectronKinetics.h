//
// Created by daan on 13-5-19.
//

#ifndef LOKI_CPP_ELECTRONKINETICS_H
#define LOKI_CPP_ELECTRONKINETICS_H

#include "Setup.h"
#include "EedfGasMixture.h"
#include "Grid.h"
#include "WorkingConditions.h"
#include "LinearAlgebra.h"

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

        Vector g_c;

        Matrix elasticMatrix,
                continuousMatrix;

        Vector eedf;

        void plot(const std::string &title, const std::string &xlabel, const std::string &ylabel,
                  const Vector &x, const Vector &y);

    public:
        explicit ElectronKinetics(const ElectronKineticsSetup &setup, const WorkingConditions *workingConditions);

        ~ElectronKinetics() = default;

        // Copying this object is not allowed.
        ElectronKinetics(const ElectronKinetics &other) = delete;
    };
}


#endif //LOKI_CPP_ELECTRONKINETICS_H
