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

// TODO: comment ElectronKinetics class

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

        Matrix elasticMatrix,
                fieldMatrix,
                CARMatrix,
                continuousMatrix,
                inelasticMatrix;

        Vector g_c, g_E, g_CAR;

        Vector eedf;

    public:
        explicit ElectronKinetics(const ElectronKineticsSetup &setup, const WorkingConditions *workingConditions);

        ~ElectronKinetics() = default;

        // Copying this object is not allowed.
        ElectronKinetics(const ElectronKinetics &other) = delete;

        void solve();

    private:

        void evaluateMatrix();

        void evaluateElasticOperator();

        void evaluateFieldOperator();

        void evaluateCAROperator();

        void evaluateInelasticOperators();

        void plot(const std::string &title, const std::string &xlabel, const std::string &ylabel,
                  const Vector &x, const Vector &y);
    };
}


#endif //LOKI_CPP_ELECTRONKINETICS_H
