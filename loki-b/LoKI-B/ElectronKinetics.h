//
// Created by daan on 13-5-19.
//

#ifndef LOKI_CPP_ELECTRONKINETICS_H
#define LOKI_CPP_ELECTRONKINETICS_H

#include "LoKI-B/EedfCollision.h"
#include "LoKI-B/EedfGasMixture.h"
#include "LoKI-B/Event.h"
#include "LoKI-B/Grid.h"
#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/MacroscopicQuantities.h"
#include "LoKI-B/Power.h"
#include "LoKI-B/Setup.h"
#include "LoKI-B/WorkingConditions.h"

// TODO: comment ElectronKinetics class

namespace loki
{

using ResultEvent =
    Event<const Grid, const Vector, const WorkingConditions, const Power, const std::vector<EedfGas *>,
          const SwarmParameters, const std::vector<RateCoefficient>, const std::vector<RateCoefficient>, const Vector>;

class ElectronKinetics
{
  public:
    ElectronKinetics(const ElectronKineticsSetup &setup, WorkingConditions *workingConditions);
    ElectronKinetics(const json_type &cnf, WorkingConditions *workingConditions);
    // Copying this object is not allowed.
    ElectronKinetics(const ElectronKinetics &other) = delete;
    // use the detaul destructor
    ~ElectronKinetics() = default;

    ResultEvent obtainedNewEedf;

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

    WorkingConditions *workingConditions;
    Grid grid;
    EedfGasMixture mixture;

    EedfType eedfType;
    /** \todo This is not used anywhere at the moment.
     *  This is used only when EEDFType is equal to 'prescribed'.
     *  This is not yet implemented in the C++ version, but present
     *  in the MATLAB version.
     */
    uint8_t shapeParameter;
    double mixingParameter;
    double maxEedfRelError;
    double maxPowerBalanceRelError;
    IonizationOperatorType ionizationOperatorType;
    GrowthModelType growthModelType;

    double CIEff{0.};
    double alphaRedEff{0.};

    // support for elastic contributions
    SparseMatrix elasticMatrix;
    Vector g_c;

    /** \todo Clarify: does inelastic mean particle-conserving inelastic only
     *                 (no ionization, attachment)? So only excitation? If so,
     *                 only electronic, or also vibrational, rotational?
     */
    Matrix inelasticMatrix;

    // Support for ionization.
    /** \todo Document whre/when ionConservativeMatrix is used. It is also used
     *        to obtain an initial guess if iterations are done, it seems.
     */
    Matrix ionConservativeMatrix;
    /** This appears to be NOT used for IonizationOperatorType::conservative
     */
    Matrix ionizationMatrix;

    // Support for electron attachment
    /** \todo Always used when attachment processes are present.
     *  Unlike ionConservativeMatrix, this does not depend on a setting, it seems.
     *  is that intended?
     */
    Matrix attachmentConservativeMatrix;
    /** \todo Always used when attachment processes are present.
     *  Unlike ionizationMatrix, this does not depend on a setting, it seems.
     *  is that intended?
     */
    SparseMatrix attachmentMatrix;

    /* Support for the E-field terms.
     * Depending on the growth model, more fields are used to handle these terms
     * (fieldMatrixSpatGrowth, g_fieldSpatialGrowth or fieldMatrixTempGrowth, g_fieldTemporalGrowth).
     */
    SparseMatrix fieldMatrix;
    Vector g_E;

    /// CARMatrix, g_CAR are relevant only when CAR gases are present (same for power.carXXX)
    SparseMatrix CARMatrix;
    Vector g_CAR;

    // variables related to spatial growth
    SparseMatrix fieldMatrixSpatGrowth;
    SparseMatrix ionSpatialGrowthD;
    SparseMatrix ionSpatialGrowthU;
    Vector g_fieldSpatialGrowth;

    // variables related to temporal growth
    SparseMatrix fieldMatrixTempGrowth;
    SparseMatrix ionTemporalGrowth;
    Vector g_fieldTemporalGrowth;

    /// \todo alphaEE, BAee, A,B are relevant only when EE collisions are configured
    bool includeEECollisions;
    /// \todo See if/when the 0-initialization is needed
    double alphaEE{0.};
    Matrix BAee;
    Vector A;
    Vector B;

    Matrix boltzmannMatrix;
    // the EEDF
    Vector eedf;
    /** \todo Note that this is essentially output only.
     *  This is evaluated at the end of solve() by a call to
     *  evaluateFirstAnisotropy(), then passed on to
     *  obtainedNewEedf.emit(), so output can be generated.
     *  This is just one example of variables of this type.
     *  It would be good to separate the variables that are
     *  needed for the calculation from those that are more
     *  of the postprecessing type.
     *  Another example is swarmParameters. Also power is mostly
     *  like this, with the exception of power.electronElectron
     *  and power.reference, which are used as part of a convergence
     *  criterium in solveEEColl().
     */
    Vector firstAnisotropy;
    Power power;
    SwarmParameters swarmParameters;

    /** \todo superElasticThresholds is only used in the disabled code path
     *  in invertMatrix that uses LinAlg::hessenbergReductionPartialPiv.
     *  It seems that all this should be removed or at least disabled as long
     *  as that code is not active.
     */
    std::vector<uint32_t> superElasticThresholds;

    /* NOTE: the following three are not configuration parameters, but the
     * result of introspection of the reaction lists. The results also depend
     * on uMax.
     * \bug It seems that the values are not always reset to false before making
     *      them conditionally true during introspection. (That is needed when uMax
     *      changes.)
     */
    bool includeNonConservativeIonization{false};
    bool includeNonConservativeAttachment{false};
    bool hasSuperelastics{false};
};

} // namespace loki

#endif // LOKI_CPP_ELECTRONKINETICS_H
