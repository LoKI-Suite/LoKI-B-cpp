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
#include "Power.h"

// TODO: comment ElectronKinetics class

namespace loki
{
using namespace Enumeration;

class ElectronKinetics
{
    EedfType eedfType;
    uint8_t shapeParameter;
    double mixingParameter;
    double maxEedfRelError;
    double maxPowerBalanceRelError;
    IonizationOperatorType ionizationOperatorType;
    GrowthModelType growthModelType;
    bool includeEECollisions{false},
        includeNonConservativeIonization{false},
        includeNonConservativeAttachment{false},
        hasSuperelastics{false};

    double CIEff{0.}, alphaRedEff{0.};

    const WorkingConditions *workingConditions;

    Grid grid;

    EedfGasMixture mixture;

    Matrix elasticMatrix,
        fieldMatrix,
        CARMatrix,
        continuousMatrix,
        inelasticMatrix,
        ionConservativeMatrix,
        ionizationMatrix,
        ionSpatialGrowthD,
        ionSpatialGrowthU,
        attachmentMatrix,
        attachmentConservativeMatrix,
        fieldMatrixSpatGrowth,
        fieldMatrixTempGrowth,
        ionTemporalGrowth,
        Aee,
        Bee;

    Vector g_c, g_E, g_CAR, g_fieldSpatialGrowth, g_fieldTemporalGrowth;

    Vector eedf;

    Power power;

    std::vector<uint32_t> superElasticThresholds;

public:
    explicit ElectronKinetics(const ElectronKineticsSetup &setup, const WorkingConditions *workingConditions);

    ~ElectronKinetics() = default;

    // Copying this object is not allowed.
    ElectronKinetics(const ElectronKinetics &other) = delete;

    void solve();

private:
    void evaluateMatrix();

    void invertLinearMatrix();

    void invertMatrix(Matrix &matrix);

    void evaluateElasticOperator();

    void evaluateFieldOperator();

    void evaluateCAROperator();

    void evaluateInelasticOperators();

    void evaluateIonizationOperator();

    void evaluateAttachmentOperator();

    void mixingDirectSolutions();

    void solveSpatialGrowthMatrix();

    void solveTemporalGrowthMatrix();

    void solveEEColl();

    void evaluatePower(bool isFinalSolution);

    void plot(const std::string &title, const std::string &xlabel, const std::string &ylabel,
              const Vector &x, const Vector &y);
};
} // namespace loki

#endif //LOKI_CPP_ELECTRONKINETICS_H
