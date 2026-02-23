/** \file
 *
 *  Unit tests for the ConvectionDiffusionOperator
 *
 *  \author Jan van Dijk
 *  \date   February 2026
 */

#include "LoKI-B/Grid.h"
#include "LoKI-B/EedfUtilities.h"
#include "LoKI-B/Operators.h"
#include <cassert>
#include <cmath>
#include <fstream>
#include <ios>

#include "tests/TestUtilities.h"

namespace loki {

/** This convective-diffusive operator manages a std::vector of pointers to
 *  external convective-diffusive operators. It calculates the base class'
 *  convection and diffusive coefficients as the sums of the contributions.
 *  These sums are re-calculated when member update() is called. Terms can
 *  be added with member register_term.
 *
 *  \author Jan van Dijk
 *  \date February 2026
 */
class ConvectionDiffusionTerms : public ConvectionDiffusionOperator
{
public:
    using Terms = std::vector<const ConvectionDiffusionOperator*>;
    using size_type = Terms::size_type;
    const Terms& terms() const { return m_terms; }
    void register_term(ConvectionDiffusionOperator& op)
    {
        m_terms.push_back(&op);
    }
    void update()
    {
        if (m_terms.empty())
        {
            this->m_conv_coeff.resize(0);
            this->m_diff_coeff.resize(0);
            return;
        }
        auto c = m_terms.begin();
        this->m_conv_coeff = (*c)->conv_coeff();
        this->m_diff_coeff = (*c)->diff_coeff();
        assert(this->m_conv_coeff.size()==this->m_diff_coeff.size());
        for ( ++c ; c != m_terms.end(); ++c)
        {
            assert((*c)->conv_coeff().size()==this->m_conv_coeff.size());
            assert((*c)->diff_coeff().size()==this->m_diff_coeff.size());
            this->m_conv_coeff += (*c)->conv_coeff();
            this->m_diff_coeff += (*c)->diff_coeff();
        }
    }
private:
    Terms m_terms;
};

/** This namespace contains various types for discretizing energy-space
 *  fluxes that are represented by ConvectionDiffusionOperator objects, and
 *  their energy-derivatives.
 *
 *  Grid layout:
 *
 *  x  |  x  |  x
 *  W  w  P  e  E
 *
 *  For f = w,e, with B(w)=W, A(w)=P, B(e)=P, A(e)=E:
 *
 *  f_f = (1-c_f)*f_B + c_f*f_A, c_f = (u_f-u_B)/(u_A-u_B).
 *
 *  Gamma_f = C_f*f_f - D_f*[df/du]_f
 *          = C_f*((1-c_f)*f_B + c_f*f_A) - (D_f/du_f)*(f_A-f_B)
 *          = (C_f*(1-c_f) + D_f/du_f)*f_B - (-C_f*c_f + D_f/du_f)*f_A
 *         := B_f*f_B - A_f*f_A,
 *
 *  with
 *
 *         B_f = C_f*(1-c_f) + D_f/du_f, A_f = -C_f*c_f + D_f/du_f.
 *
 *  Gamma_e-Gamma_w = (B_e*f_P - A_e*f_E) - (B_w*f_W - A_w*f_P)
 *                  = -B_wf_W + A_w*f_P + B_e*f_P - A_e*f_E.
 *
 *  At the eastern interface, we can introduce a ghost point f_E at
 *  position u_max. The flux is then given by
 *
 *         Gamma_bnd = B_bnd*f_P - A_bnd*f_E,
 *
 *  and setting this to 0 results in f_E = f_P*B_f/A_f. Here f_P is
 *  the value in the last real cell. With this choice for the ghost
 *  point energy, c_f=1 and du_f = du_Pe = du_we/2.
 *
 *  Below we discretize -(Gamma_e-Gamma_w)/du_we --- note the extra
 *  minus sign.
 *
 *  \todo Note that this renders the matrix NEGATIVE semi-definite.
 *
 *  \todo Check the upper boundary flux discretization (ghost point location,
 *  in particular, see if that is used consistently).
 *
 *  \author Jan van Dijk
 *  \date February 2026
 */
namespace Schemes
{
    /** The meaning of these coefficients for some face index f is that
     *  Gamma_f = B*f_B - A*f_A, where f_B and f_A are the field values
     *  before and after the face.
     */
    struct LocalFluxCoefficients
    {
        double A;
        double B;
    };
    struct CD
    {
        /** For the central difference scheme, the partial flux does not depend
         *  on the total flux. This means that the resulting coefficients do not
         *  depend on \a sum, except at the upper boundary, where the field value
         *  follows from the requirement that the *total* flux is zero.
         */
        static LocalFluxCoefficients calc_coefs(
            const Grid& grid,
            const ConvectionDiffusionOperator& term,
            const ConvectionDiffusionOperator& sum,
            Grid::Index k)
        {
            const double du_Lf = grid.duCell(k-1) / 2;
            const double du_LH = grid.duNode(k);
            const double cf = du_Lf / du_LH;
            const double A = term.diff_coeff()[k]/du_LH - term.conv_coeff()[k]*cf;
            // todo: (?) rewrite as const double B = A + term.conv_coeff()[k];
            const double B = term.diff_coeff()[k]/du_LH + term.conv_coeff()[k]*(1-cf);
            if (k == grid.getNodes().size() - 1)
            {
                /* At the upper boundary, the total flux (represented by 'sum')
                 * is given by Bsum*f_B - Asum*f_g = 0. Then
                 * Gamma = B*f_B - A*f_g = B*f_B - A*f_B*(Bsum/Asum) = (B - A*(Bsum/Asum))*f_B.
                 */
                const double sum_A = sum.diff_coeff()[k]/du_LH - sum.conv_coeff()[k]*cf;
                // todo (?) rewrite as const double sum_B = sum_A + sum.conv_coeff()[k];
                const double sum_B = sum.diff_coeff()[k]/du_LH + sum.conv_coeff()[k]*(1-cf);
                return { std::numeric_limits<double>::quiet_NaN(), B - A*sum_B/sum_A };
            }
            return { A, B };
        }
    };
    struct SG
    {
        static constexpr double bernoulli(double x)
        {
            return std::abs(x) < 1e-9 ? 1-x/2 : x/std::expm1(x);
        }
        /** \a k is a face index.
         */
        static LocalFluxCoefficients calc_coefs(
            const Grid& grid,
            const ConvectionDiffusionOperator& term,
            const ConvectionDiffusionOperator& sum,
            Grid::Index k)
        {
            const double du_LH = grid.duNode(k);
            const double Pe = sum.conv_coeff()[k]*du_LH/sum.diff_coeff()[k];
            const double ber = bernoulli(Pe);
            if (&term==&sum)
            {
                /* This means that we are discretizing the total flux. This
                 * allows some simplifications of the evaluation, but the result
                 * is the same as (barring round-off errors).
                 */
                if (k == grid.getNodes().size() - 1)
                {
                    return { std::numeric_limits<double>::quiet_NaN(), 0.0 };
                }
                const double A = term.diff_coeff()[k]/du_LH*ber;
                const double B = term.diff_coeff()[k]/du_LH*(ber+Pe);
                return { A, B };
            }
            // Handle the general case: term is a partial flux.
            const double du_Lf = grid.duCell(k-1) / 2;
            const double cf = du_Lf / du_LH;
            /* NOTE: for small P we use the approximation
             * expm1(Pc)/expm1(P) = (1 + Pc + ... -1)/(1 + P + ... -1) = c + ...
             */
            const double expm1_ratio = std::abs(Pe)<1e-9 ? cf : std::expm1(Pe*cf)/std::expm1(Pe);
            const double A = term.diff_coeff()[k]/du_LH*ber*std::exp(Pe*cf) - term.conv_coeff()[k]*expm1_ratio;
            const double B = A + term.conv_coeff()[k];
            if (k == grid.getNodes().size() - 1)
            {
                /* At the upper boundary, the total flux (represented by 'sum')
                 * is given by Bsum*f_B - Asum*f_g = 0. Then
                 * Gamma = B*f_B - A*f_g = B*f_B - A*f_B*(Bsum/Asum) = (B - A*(Bsum/Asum))*f_B.
                 */
                const double sum_A = sum.diff_coeff()[k]/du_LH*ber;
                const double sum_B = sum.diff_coeff()[k]/du_LH*(ber+Pe);
                return { std::numeric_limits<double>::quiet_NaN(), B - A*sum_B/sum_A };
            }
            return { A, B };
        }
    };

} // namespace Schemes

/** Flux contribution, using the SchemeTraits.
 *  This uses the contribution's convection and diffusion coefficients,
 *  and the Peclet number that is based on the sums of the convective and
 *  diffusive contributions.
 */
template <class SchemeTraits>
void discretize_dflux_du(Matrix& mat, const Grid& grid, const ConvectionDiffusionOperator& term, const ConvectionDiffusionOperator& sum)
{
    mat.fill(0.0);
    /** \todo We calculate every face twice. This can be fixed by changing
     *  this in a loop over the faces, and handle the contributions to the
     *  equations in both the lower and upper adjacent cells.
     */
    for (Grid::Index k = 1; k < grid.nCells(); ++k)
    {
        mat.coeffRef(k, k) = 0;

        const double du_we = grid.duCell(k);
        // The flux is 0 at the lower boundary.
        if (k > 0)
        {
            const auto coefs = SchemeTraits::calc_coefs(grid,term,sum,k);
            mat.coeffRef(k, k - 1)  = +coefs.B / du_we;
            mat.coeffRef(k, k)     += -coefs.A / du_we;
        }
#define LOKI_DISCRETIZE_BOUNDARY_FLUX 1
#if LOKI_DISCRETIZE_BOUNDARY_FLUX
        /* In general, the flux is given by Gamma = B*f_B - A*f_A. At the
         * upper boundary, Gamma = B*f_B instead. (calc_coefs modifies
         * coefs.B by eliminating value f_A in the ghost point outide of the grid).
         */
        const auto coefs = SchemeTraits::calc_coefs(grid,term,sum,k+1);
        mat.coeffRef(k, k)     += -coefs.B / du_we;
        if (k<grid.nCells()-1)
        {
            mat.coeffRef(k, k + 1)  = +coefs.A / du_we;
        }
#else
        if (k < grid.nCells() - 1)
        {
            const auto coefs = SchemeTraits::calc_coefs(grid,term,sum,k+1);
            mat.coeffRef(k, k)     += -coefs.B / du_we;
            mat.coeffRef(k, k + 1)  = +coefs.A / du_we;
        }
#endif
    }
}

/** Discretize the total flux, using the SchemeTraits.
 */
template <class SchemeTraits>
void discretize_dflux_du(Matrix& mat, const Grid& grid, const ConvectionDiffusionOperator& term)
{
    discretize_dflux_du<SchemeTraits>(mat,grid,term,term);
}

template <class SchemeTraits>
void evaluate_flux(Vector& flux, const Grid& grid, const Vector& eedf, const ConvectionDiffusionOperator& term, const ConvectionDiffusionOperator& sum)
{
    flux[0] = 0.0;
    for (Grid::Index k = 1; k < grid.getNodes().size()-1; ++k)
    {
            const auto coefs = SchemeTraits::calc_coefs(grid,term,sum,k);
            flux[k] = coefs.B*eedf[k-1] - coefs.A*eedf[k];
    }
    const Grid::Index k = grid.getNodes().size()-1;
#if LOKI_DISCRETIZE_BOUNDARY_FLUX
    const auto coefs = SchemeTraits::calc_coefs(grid,term,sum,k);
    flux[k] = coefs.B*eedf[k-1];
    assert(coefs.A==std::numeric_limits<double>::quiet_NaN());
#else
    flux[k] = 0.0;
#endif
}

} // namespace loki

/* Call solveEEDF to solve the 'Boltzmann matrix' mat, calculate the mean
 * energy of the resulting EEDF and calculate the Druyvesteyn distribution
 * for this mean energy (this is the analytical solution in case Tg=0).
 * Write the energies, eedf and Druyvesteyn eedfs to file f_<suffix>.dat
 */
loki::Vector solveCase(const loki::Grid& grid, loki::Matrix& mat, const std::string& suffix)
{
    std::cout << "Running case '" << suffix << "'." << std::endl;
    std::cout << " * Mold: " << std::endl << mat << std::endl;
    loki::Vector eedf(grid.nCells());
    loki::solveEEDF(eedf,mat,grid);

    const double eps = loki::getMeanEnergy(eedf,grid);
    const double T_eV = eps*2.0/3.0;
    // 2: Druyvesteyn
    const loki::Vector eedf_an = loki::makePrescribedEDF(grid, 2, T_eV);

    std::cout << " * T = " << T_eV << " eV" << std::endl;

    std::ofstream ofs("f_" + suffix + ".dat");
    for (loki::Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        ofs << grid.getCells()[k] << '\t' << eedf[k] << '\t' << eedf_an[k] << std::endl;
    }
    return eedf;
}

template <class SchemeT>
void test_scheme(const loki::Grid& grid, const loki::ConvectionDiffusionTerms& conv_diff_terms, const std::string& scheme_name)
{
    using namespace loki;

    Matrix mat(grid.nCells(),grid.nCells());

    /* 1. Do the discretization of the total flux.
     */

    mat.fill(0.0);
    discretize_dflux_du<SchemeT>(mat,grid,conv_diff_terms);
    solveCase(grid,mat,scheme_name);

    /* 2. Do the discretization of the individual flux contributions, store the
     * matrices. Aggregate the results afterwards to discretize the total flux,
     * then evaluate the individual power terms.
     */
    {

    std::vector<Matrix> mats(conv_diff_terms.terms().size());
    Matrix mat_term(grid.nCells(),grid.nCells());
    for (typename ConvectionDiffusionTerms::size_type t=0; t!=conv_diff_terms.terms().size(); ++t)
    {
        mats[t].resize(grid.nCells(),grid.nCells());
        mats[t].fill(0.0);
        discretize_dflux_du<SchemeT>(mats[t],grid,*conv_diff_terms.terms()[t],conv_diff_terms);
        mat += mats[t];
    }
    const auto eedf = solveCase(grid,mat,scheme_name+"_sep");

    std::vector<Vector> fluxes(conv_diff_terms.terms().size());
    Vector flux_sum(grid.getNodes().size());
    flux_sum.fill(0.0);
    for (typename ConvectionDiffusionTerms::size_type t=0; t!=conv_diff_terms.terms().size(); ++t)
    {
        fluxes[t].resize(grid.getNodes().size());
        evaluate_flux<SchemeT>(fluxes[t],grid,eedf,*conv_diff_terms.terms()[t],conv_diff_terms);
        flux_sum += fluxes[t];
    }
    Vector flux_total(grid.getNodes().size());
    evaluate_flux<SchemeT>(flux_total,grid,eedf,conv_diff_terms,conv_diff_terms);
    std::ofstream ofs("flux_"+scheme_name+"_sep.dat");
    for (loki::Grid::Index k = 0; k < grid.getNodes().size(); ++k)
    {
        ofs << grid.getNodes()[k];
        for (const auto& f : fluxes)
        {
            ofs << '\t' << f[k];
        }
        ofs << '\t' << flux_sum[k] << '\t' << flux_total[k] << std::endl;
    }

    /* NOTE: at this point, the first rows of the matrices are modified to
     * incorporate the normalization constraint. But we do not care, since
     * the row is multiplied with a zero energy value in the calculation of
     * the power terms using the expressions below.
     *
     * This statement of the power terms is explained in the cpp_notes document.
     */
    const Vector gamma_u_du = SI::gamma*grid.duCells().array()*grid.getCells().array();
    std::cout << "Power balance --- volumetric power/(n_e*N) in eV*m^3/s:" << std::endl;
    double P_sum=0;
    for (typename ConvectionDiffusionTerms::size_type t=0; t!=conv_diff_terms.terms().size(); ++t)
    {
        const double P_term = gamma_u_du.dot(mats[t]*eedf);
        std::cout << "Term #" << t << ": " << std::showpos << P_term << std::endl;
        P_sum += P_term;
    }
    std::cout << "    Sum: " << std::showpos << P_sum << std::endl;

    }
}

int main()
{
    using namespace loki;
    //const double eps = std::numeric_limits<double>::epsilon();

    constexpr Grid::Index Nc = 1001;
    constexpr double u_max = 1.4; // eV
    const Grid grid(Nc,u_max);

    const double Tg = 0; // K
    const double WoN = 0.0;
    const double CIEff = 0.0;
    const double EoN = 10e-21; // Vm^2

    Vector sigma_f(grid.getNodes().size());
    sigma_f.fill(1e-19);

    ElasticOperator op_elast;
    FieldOperator op_field(grid);

    {
        /* 1. The traditional solution method: call evaluate on the operators
         * to discretize the matrix contributions.
         */
        SparseMatrix mat_elast(grid.nCells(),grid.nCells());
        SparseMatrix mat_field(grid.nCells(),grid.nCells());

        op_elast.evaluate(grid,sigma_f,Tg,mat_elast);
        op_field.evaluate(grid,sigma_f,EoN,WoN,CIEff,mat_field);

        Matrix mat = mat_elast + mat_field;

        const auto eedf = solveCase(grid,mat,"old");

        double P_elast, P_field;
        double dummy;
        op_elast.evaluatePower(grid,eedf,Tg,P_elast,dummy,dummy);
        /* Here we pass E/N=1. This is a result of the awkward scaling of the
         * field term: when evaluate is called with 1, we need to apply the
         * factor (E/N)^2 here; when that is called with the actual E/N value,
         * we pass 1 here to prevent a second scaling with (E/N)^2.
         */
        op_field.evaluatePower(grid,eedf,1.0,P_field);
        std::cout << "Power balance --- volumetric power/(n_e*N) in eV*m^3/s:" << std::endl;
        std::cout << " Elastic: " << std::showpos << P_elast << std::endl;
        std::cout << " Field:   " << std::showpos << P_field << std::endl;
        std::cout << " Net:     " << std::showpos << (P_elast+P_field) << std::endl;

    }

    /* the following calculations solve the same problem, but using the
     * ConvectionDiffusionDiscretizer class. First we register the (same)
     * elastic and field operators, next we call updateCD on both to update
     * their convection and diffusion coefficients. Then we call update()
     * on the discretizer object, which will update the summed convection
     * and diffusion coefficients and the grid Peclet number.
     */
    ConvectionDiffusionTerms conv_diff_terms;
    conv_diff_terms.register_term(op_elast);
    conv_diff_terms.register_term(op_field);

    op_elast.updateCD(grid,sigma_f,Tg);
    op_field.updateCD(grid,sigma_f,EoN,WoN,CIEff);

    conv_diff_terms.update();

    test_scheme<Schemes::CD>(grid, conv_diff_terms, "cd");
    test_scheme<Schemes::SG>(grid, conv_diff_terms, "sg");

    return 0;
}
