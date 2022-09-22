#ifndef LOKI_CPP_OPERATORS_H
#define LOKI_CPP_OPERATORS_H

#include "LoKI-B/Gas.h"
#include "LoKI-B/Grid.h"
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


} // namespace loki

#endif // LOKI_CPP_OPERATORS_H
