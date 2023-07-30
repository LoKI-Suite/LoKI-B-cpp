#ifndef LOKI_CPP_OPERATORS_H
#define LOKI_CPP_OPERATORS_H

#include "LoKI-B/Gas.h"
#include "LoKI-B/Grid.h"
#include "LoKI-B/EedfMixture.h"
#include <vector>

namespace loki {

    /** When comparing this with Tejero2019, realize that in that paper,
     * equation 6c, the following symbols are used for a gas k:
     *
     * B_k: rotationalConstant
     * Q_{k,au}: quadruple moment in (atomic) units e*a_0^2. Here e is
     *   the elementary charge and a_0 the Bohr radius. (NOTE that the
     *   variable electricQuadrupoleMoment in the code is in SI units Cm^2.)
     * sigma_{0,k} = (8./15)*pi*Q_{k,au}^2*a_0^2, see \cite Tejero below
     * equation 6d, \cite Ridenti below equation 8b or Gerjuoy and Stein,
     * equation 20.
     *
     * For mixtures, the terms B_k*sigma_k in the expression for g must be
     * weighted with the molar fractions. The code first calculates this weighted
     * sum sigma0B, which gives
     *
     *   sum_k chi_k*B_k*sigma_k = (8./15)*pi*a_0^2* sum_k chi_k*B_k*Q_{k,au}^2
     *
     * NOTE: in the first public release of LoKI-B (MATLAB), Q_{k,au} was used
     *       instead of Q_{k,au}^2. That has been fixed in version 2.0.0.
     *
     *  \todo In the code all the G's are divided by N*sqrt(2*e/m_e),
     *  compared to the LoKI-B paper \cite Tejero2019, it seems.
     *  That explains why, in the code below, you see gas->fraction,
     *  whereas in the paper you see N_k. It would be good to have a
     *  document where the equations are written *exactly* as in the code.
     *  Also g is defined without the minus sign that appears in the definition
     *  in the paper. All in all, CARmatrix seems to be defined such that
     *  [CARmatrix]*[f] is an approximation of -(1/(N*sqrt(2*e/m_e))dG_CAR/du.
     */
    class CAROperator
    {
    public:
        using CARGases = std::vector<const Gas*>;
        CAROperator(const CARGases& cg);
        /// updates member g
        void evaluate(const Grid& grid);
        /// updates member g, then the CAR matrix \a mat
        void evaluate(const Grid& grid, double Tg, SparseMatrix& mat);
        void evaluatePower(const Grid& grid, const Vector& eedf, double Tg, double& net, double& gain, double& loss) const;
        Vector g;
        const CARGases carGases;
    };

    class ElasticOperator
    {
    public:
        ElasticOperator();
        /// updates member g
        void evaluate(const Grid& grid, const Vector& elasticCrossSection);
        /// updates member g, then the elastic matrix \a mat
	void evaluate(const Grid& grid, const Vector& elasticCrossSection, double Tg, SparseMatrix& mat);
        void evaluatePower(const Grid& grid, const Vector& eedf, double Tg, double& net, double& gain, double& loss) const;
        Vector g;
    };

    class FieldOperator
    {
    public:
        FieldOperator(const Grid& grid);
        /// updates member g
        void evaluate(const Grid& grid, const Vector& totalCS, double EoN, double WoN);
        /// updates member g, then the field matrix \a mat
	void evaluate(const Grid& grid, const Vector& totalCS, double EoN, double WoN, SparseMatrix& mat);
        /** \todo This expression returns the power that is absorbed from the field
         *  in the case that no temporal or spatial growth terms are present. We have
         *  to see where/when to handle the various growth scenarios.
         */
        void evaluatePower(const Grid& grid, const Vector& eedf, double& power) const;
        Vector g;
    };

    class InelasticOperator
    {
    public:
	InelasticOperator(const Grid& grid);

        /** \todo Clarify: does inelastic mean particle-conserving inelastic only
         *                 (no ionization, attachment)? So only excitation? If so,
         *                 only electronic, or also vibrational, rotational?
         */
        void evaluateInelasticOperators(const Grid& grid, const EedfMixture& mixture);
        Matrix inelasticMatrix;
        /* NOTE: the following is not a configuration parameter, but the
         * result of introspection of the reaction lists. The results also depend
         * on uMax.
         * \bug It seems that the value is not always reset to false before making
         *      it conditionally true during introspection. (That is needed when uMax
         *      changes.)
         */
        bool hasSuperelastics;
    };

    class ElectronElectronOperator
    {
    public:
        ElectronElectronOperator();
        void initialize(const Grid& grid);
	void updateABMatrices(const Grid& grid);
	void update_g_ee(const Grid& grid, const Vector& eedf, double ne, double n0);
	/** updateAB requires that updateABMatrices and update_g_ee have been called
	 *  to set up BAee and g_ee.
	 *  \todo Should update_g_ee not always be called before updateAB is called,
	 *  since Te may have changed. At present this is not the case. Check the
	 *  users of this code.
	 */
	void updateAB(const Grid& grid, const Vector& eedf);
	/** Adds the coefficients coming from the discretization of the term
         *  -(1/N*gamma)*dG_ee/du to matrix \a M.
	 */
        void discretizeTerm(Matrix& M, const Grid& grid) const;
        void evaluatePower(const Grid& grid, const Vector& eedf, double& power) const;
        double g_ee;
        Matrix BAee;
        Vector A;
        Vector B;
    };

    class IonizationOperator
    {
    public:
        IonizationOperator(IonizationOperatorType type);
        void evaluateIonizationOperator(const Grid& grid, const EedfMixture& mixture);
        const IonizationOperatorType ionizationOperatorType;
        /** \todo Document whre/when ionConservativeMatrix is used. It is also used
         *        to obtain an initial guess if iterations are done, it seems.
         */
        Matrix ionConservativeMatrix;
        /** This appears to be NOT used for IonizationOperatorType::conservative
         */
        Matrix ionizationMatrix;
        /* NOTE: the following is not a configuration parameter, but the
         * result of introspection of the reaction lists. The results also depend
         * on uMax.
         * \bug It seems that the value is not always reset to false before making
         *      it conditionally true during introspection. (That is needed when uMax
         *      changes.)
         */
        bool includeNonConservativeIonization;
    };

    class AttachmentOperator
    {
    public:
        // Support for electron attachment
        void evaluateAttachmentOperator(const Grid& grid, const EedfMixture& mixture);

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
        /* NOTE: the following is not a configuration parameter, but the
         * result of introspection of the reaction lists. The results also depend
         * on uMax.
         * \bug It seems that the value is not always reset to false before making
         *      it conditionally true during introspection. (That is needed when uMax
         *      changes.)
         */
        bool includeNonConservativeAttachment{false};
    };

} // namespace loki

#endif // LOKI_CPP_OPERATORS_H
