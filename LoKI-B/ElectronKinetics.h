//
// Created by daan on 13-5-19.
//

#ifndef LOKI_CPP_ELECTRONKINETICS_H
#define LOKI_CPP_ELECTRONKINETICS_H

#include "LoKI-B/EedfCollisions.h"
#include "LoKI-B/EedfMixture.h"
#include "LoKI-B/Event.h"
#include "LoKI-B/Grid.h"
#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/MacroscopicQuantities.h"
#include "LoKI-B/Power.h"
#include "LoKI-B/WorkingConditions.h"
#include "LoKI-B/Operators.h"

// TODO: comment ElectronKinetics class

namespace loki
{

using ResultEvent =
    Event<const Grid&, const Vector&, const WorkingConditions&, const Power&, const EedfCollisionDataMixture&,
          const SwarmParameters&, const Vector*>;

class ElectronKinetics
{
protected:
    ElectronKinetics(const std::filesystem::path &basePath, const json_type &cnf, WorkingConditions *workingConditions);
    // Copying this object is not allowed.
    ElectronKinetics(const ElectronKinetics &other) = delete;
public:
    // use the detaul destructor
    virtual ~ElectronKinetics() = default;
    /** \todo Implement the 'rule of big 5': also delete the move constructor and
     *  move and copy assignment.
     */
    ResultEvent obtainedNewEedf;

    /** solve the Boltzmann equation.
     *  After completion, the power terms, rate coefficients, swarm parameters
     *  and the first anisotropic term f1 are evaluated and the obtainedNewEedf
     *  event is fired.
     */
    void solve();

    const Grid &grid() const { return m_grid; }
protected:
    /** Calls updateMaxEnergy(uMax) on the grid. This function allows us to
     *  keep the grid member private and expose only a non-constant reference
     *  to the grid (via member grid()).
     */
    void updateMaxEnergy(double uMax);
    void updateMaxEnergyNonuniform(double uMax);
    virtual void doSolve()=0;

    /// \todo See what can be made private. Introduce accessors where necessary.

    /** Given two numbers, calculate how many decades |v2| is smaller than |v1|.
     *  The result is calculated as log10(|v1/v2|). Note that calcDecades(0,0)=NaN.
     *  Furthermore, calcDecades(v1,v2)=-calcDecades(v2,v1) and in particular,
     *  for non-zero v we have calcDecades(0,v)=-Inf and calcDecades(v,0)=+Inf;
     *
     *  Loki-B uses this function to evaluate the dynamic range of the EEDF, and
     *  calculate the energy grid accordingly ('smart grid').
     *
     *  Note that the abs is needed even if v1,v2 are non-negative. In case of
     *  underflow, LoKI-B may produce v1>0 and v2=-0. Such signed zero would
     *  result in v1/v2=-Inf, for which case log(v1/v2)=NaN, where +Inf is desired.
     *  The old implementation (still present in the MATLAB code) was written
     *  as log10(v1)-log10(v2). Note that using the mathematical identity
     *  log10(a)-log10(b) = log10(a/b) could NOT be used, since for positive a
     *  and b=+/-0.0 the LHS produces log10(a)-(-Inf) = Inf, whereas the right
     *  hand side would produce log10(+/- Inf) and would depend on the sign of
     *  the zero: log10(Inf)=Inf, but log10(-Inf)=NaN.
     *  By using the abs, we can solve this problem and need to evaluate only
     *  a single logarithm.
     *
     *  \author Jan van Dijk
     *  \date   December 2020
     */
    static double calcDecades(double v1, double v2)
    {
        return std::log10(std::abs(v1/v2));
    }
    WorkingConditions *m_workingConditions;
private:
    Grid m_grid;
protected:
    EedfMixture mixture;

    ElasticOperator elasticOperator;

    /* Support for the E-field terms.
     * Depending on the growth model, more fields are used to handle these terms
     * (fieldMatrixSpatGrowth, g_fieldSpatialGrowth or fieldMatrixTempGrowth, g_fieldTemporalGrowth).
     */
    FieldOperator fieldOperator;

    InelasticOperator inelasticOperator;

    // support for CAR processes.
    std::unique_ptr<CAROperator> carOperator;

    // the EEDF
    Vector eedf;

    // storage of the power terms
    Power power;

    SwarmParameters swarmParameters;

    /** \todo superElasticThresholds is only used in the disabled code path
     *  in invertMatrix that uses LinAlg::hessenbergReductionPartialPiv.
     *  It seems that all this should be removed or at least disabled as long
     *  as that code is not active.
    std::vector<uint32_t> superElasticThresholds;
     */

private:

    /** Carry out initialization tasks. This function is called by both
     *  constructor overloads after the argument-type-specific bits have
     *  been done (those that depend on WorkingCondisions or JSON).
     */
    void initialize();
};

class ElectronKineticsBoltzmann : public ElectronKinetics
{
public:
    ElectronKineticsBoltzmann(const std::filesystem::path &basePath, const json_type &cnf, WorkingConditions *workingConditions);
protected:
    /** solve the Boltzmann equation. Calls solveSingle or solveSmartGrid,
     *  depending on wether the smart grid option is enabled.
     */
    virtual void doSolve();
private:
    /// shared constructor tasks
    void initialize();
    /** solve the Boltzmann equation, taking into account only the linear terms.
     */
    void invertLinearMatrix();
    /** solve matrix*eedf=b, with b=[0], subject to the constraint that
     *  sum_i sqrt(u_i)*f[i]*du = 1.
     *  Note: the incoming matrix is singular, the first equation (that is:
     *  matrix(0,*) annd b[0] are used to encode the normalization constraint.
     */
    void invertMatrix(Matrix &matrix);
    void solveSpatialGrowthMatrix();
    void solveTemporalGrowthMatrix();
    /** Solve the EE collision matrix while other terms are fixed.
     *  This requires that eeOperator has been set up.
     */
    void solveEEColl();
    /** Update the terms that are needed later on to compose the boltzmannMatrix
     *  that need to be calculated only once after a grid update.
     *  Note that this is called only in the constructor of this class and in
     *  reply to an grid.updatedMaxEnergy event (which is triggered among others
     *  by a smart grid iteration).
     */
    void evaluateMatrix();
    /** Solve the Boltzmann equation. The function first inverts the linear
     *  problem. If non-linear terms are available (growth terms or electron
     *  electron collisions, subsequently use the 'mixing direct solutions'
     *  algorithm to find the solution.
     */
    void obtainTimeIndependentSolution();
    /** solve the Boltzmann equation by calling obtainTimeIndependentSolution()
     *  (a time-dependent alternative, available in MATLAB, should become
     *  available later.
     */
    void solveSingle();
    void solveSmartGrid();
    /** An alternative for solveSmartGrid.
     *
     *  0. introduce uM, uP and set these to uMax, the initial value of the
     *     upper boundary.
     *  1. do a calculation and calculate 'decades'.
     *  2. If necessary, widen the interval [uM,uP].
     *     - if decades is below the interval:
     *         while decades is below the interval:
     *          - double uP, set uMax=uP, solve and recalculate 'decades'
     *     - else if decades is above the interval:
     *         while decades is above the interval:
     *          - divide uM by 2, set uMax=uM, solve and recalculate 'decades'
     *  3. While decades is not in the range, use bisection:
     *      - Set uMax = (uM+uP)/2
     *      - solve and calculate decades.
     *      - if decades is below the intervael uM <- uMax, otherwise uP <- uMax
     */
    void solveSmartGrid2();
    // storage and calculation of f1 (first anisotropic term)
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
    void evaluateFirstAnisotropy();
    Vector firstAnisotropy;
    // calculation of the power terms
    void evaluatePower();
    // storage and calculation of the swarm parameters
    void evaluateSwarmParameters();

    SparseMatrix elasticMatrix;
    void evaluateFieldOperator();
    SparseMatrix fieldMatrix;
    SparseMatrix CARMatrix;

    AttachmentOperator attachmentOperator;

    Matrix boltzmannMatrix;
    // Support for ionization.
    IonizationOperator ionizationOperator;

    const std::unique_ptr<ElectronElectronOperator> eeOperator;

    double CIEff{0.};

    // use by spatial AND temporal growth
    GrowthModelType growthModelType;
    double mixingParameter;

    // code related to spatial growth
    // term arising from Manual 2.2.0 eq. 8a with the first term of 8b:
    // (note: fieldMatrix handles the second term of 8b)
    Vector g_fieldSpatialGrowth;
    SparseMatrix fieldMatrixSpatGrowth;
    SparseMatrix ionSpatialGrowthD;
    SparseMatrix ionSpatialGrowthU;
    double alphaRedEff{0.};

    // code related to temporal growth
    Vector g_fieldTemporalGrowth;
    SparseMatrix fieldMatrixTempGrowth;
    SparseMatrix ionTemporalGrowth;
    /** Tolerance settings.
     */
    double maxEedfRelError;
    double maxPowerBalanceRelError;

};

class ElectronKineticsPrescribed : public ElectronKinetics
{
public:
    ElectronKineticsPrescribed(const std::filesystem::path &basePath,const json_type &cnf, WorkingConditions *workingConditions);
protected:
    virtual void doSolve();
private:
    /// shared constructor tasks
    void initialize();
    void evaluateFieldOperator();
    void evaluateMatrix();
    // calculation of the power terms
    void evaluatePower();
    // storage and calculation of the swarm parameters
    void evaluateSwarmParameters();
    /** The parameter that controls the shape of the eedf (1 Maxwellian,
     *  2 Druyvesteyn)
     */
    const double shapeParameter;
};

} // namespace loki

#endif // LOKI_CPP_ELECTRONKINETICS_H
