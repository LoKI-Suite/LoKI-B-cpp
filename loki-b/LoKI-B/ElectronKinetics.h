//
// Created by daan on 13-5-19.
//

#ifndef LOKI_CPP_ELECTRONKINETICS_H
#define LOKI_CPP_ELECTRONKINETICS_H

#include "LoKI-B/Setup.h"
#include "LoKI-B/EedfGasMixture.h"
#include "LoKI-B/EedfCollision.h"
#include "LoKI-B/Grid.h"
#include "LoKI-B/WorkingConditions.h"
#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/Power.h"
#include "LoKI-B/MacroscopicQuantities.h"
#include "LoKI-B/Event.h"

// TODO: comment ElectronKinetics class

namespace loki
{

typedef Event<const Grid, const Vector, const WorkingConditions, const Power, const std::vector<EedfGas*>,
              const SwarmParameters, const std::vector<RateCoefficient>, const std::vector<RateCoefficient>,
              const Vector>
    ResultEvent;

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

    double CIEff{0.}, alphaRedEff{0.}, alphaEE{0.};

    WorkingConditions *workingConditions;

    Grid grid;

    EedfGasMixture mixture;

    Matrix inelasticMatrix,
        ionConservativeMatrix,
        ionizationMatrix,
        attachmentConservativeMatrix,
        BAee,
        boltzmannMatrix;

    SparseMatrix elasticMatrix,
        fieldMatrix,
        CARMatrix,
        attachmentMatrix,
        ionSpatialGrowthD,
        ionSpatialGrowthU,
        fieldMatrixSpatGrowth,
        fieldMatrixTempGrowth,
        ionTemporalGrowth;

    Vector g_c, g_E, g_CAR, g_fieldSpatialGrowth, g_fieldTemporalGrowth, A, B;

    Vector eedf, firstAnisotropy;

    Power power;

    SwarmParameters swarmParameters;

    std::vector<uint32_t> superElasticThresholds;

public:
    ElectronKinetics(const ElectronKineticsSetup &setup, WorkingConditions *workingConditions);
    ElectronKinetics(const json_type &cnf, WorkingConditions *workingConditions);

    ResultEvent obtainedNewEedf;

    ~ElectronKinetics() = default;

    // Copying this object is not allowed.
    ElectronKinetics(const ElectronKinetics &other) = delete;

    void solve();

    const Grid *getGrid();

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

    void evaluateSwarmParameters();

    void evaluateFirstAnisotropy();
};
} // namespace loki

#endif //LOKI_CPP_ELECTRONKINETICS_H
