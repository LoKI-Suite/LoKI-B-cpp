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
#include <cmath>
#include <fstream>
#include <ios>

#include "tests/TestUtilities.h"

namespace loki {

/** Grid layout:
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
 *  position u_max + du_ew[n-1]. The flux is then given by
 *
 *         Gamma_bnd = B_bnd*f_P - A_bnd*f_E,
 *
 *  and setting this to 0 results in f_E = f_P*B_f/A_f. Here f_P is
 *  the value in the last real cell. With this choice for the ghost
 *  point energy, c_f=1/2 and du_f = 2*du_Pe = du_we.
 *
 *  Below we discretize -(Gamma_e-Gamma_w)/du_we --- note the extra
 *  minus sign.
 *  \todo Note that this renders the matrix NEGATIVE semi-definite.
 */
class ConvectionDiffusionDiscretizer
{
public:
    using Contributions = std::vector<const ConvectionDiffusionOperator*>;
    ConvectionDiffusionDiscretizer(const Grid& grid)
    : m_grid(grid)
    {
    }
    const Vector& C() const { return m_C; }
    const Vector& D() const { return m_D; }
    const Vector& Pe() const { return m_Pe; }
    void update()
    {
        if (m_contribs.empty())
        {
            m_C.resize(0);
            m_D.resize(0);
            m_Pe.resize(0);
            return;
        }
        auto c = m_contribs.begin();
        m_C = (*c)->conv_coeff();
        m_D = (*c)->diff_coeff();
        for ( ++c ; c != m_contribs.end(); ++c)
        {
            m_C += (*c)->conv_coeff();
            m_D += (*c)->diff_coeff();
        }
        m_Pe = m_C.array()*m_grid.duNodes().array()/m_D.array();
    }
    void register_term(ConvectionDiffusionOperator& op)
    {
        m_contribs.push_back(&op);
    }
    /** The meaning of these coefficients for some face index f is that
     *  Gamma_f = B*f_B - A*f_A, where f_B and f_A are the field values
     *  before and after the face.
     */
    struct FluxCoefs
    {
        double A;
        double B;
    };
    /** \a k is a face index
     *  NOTE: this requires that update has been called, so the members
     *  m_C and m_D are up-to-date.
     */
    FluxCoefs calc_coefs_cd(const Vector& C, const Vector& D, Grid::Index k) const
    {
        const double du_Ww = m_grid.duCell(k-1) / 2;
        const double du_WP = m_grid.duNode(k);
        const double cf = du_Ww / du_WP;
        const double A = D[k]/du_WP - C[k]*cf;
        const double B = D[k]/du_WP + C[k]*(1-cf);
        return { A, B };
    }
    /** \todo The discretize_XXX functions in this class only differ in which
     *  calc_coefs_xxx member they call for retrieving the local coefficients
     *  A and B.
     */
    void discretize_cd(const Vector& C, const Vector& D, Matrix& mat) const
    {
        mat.fill(0.0);
        for (Grid::Index k = 1; k < m_grid.nCells(); ++k)
        {
            mat.coeffRef(k, k) = 0;

            const double du_we = m_grid.duCell(k);
            if (k > 0)
            {
                const auto coefs = calc_coefs_cd(C,D,k);
                mat.coeffRef(k, k - 1)  = +coefs.B / du_we;
                mat.coeffRef(k, k)     += -coefs.A / du_we;
            }

            if (k < m_grid.nCells() - 1)
            {
                const auto coefs = calc_coefs_cd(C,D,k+1);
                mat.coeffRef(k, k)     += -coefs.B / du_we;
                mat.coeffRef(k, k + 1)  = +coefs.A / du_we;
            }
        }
    }
    /** Discretize the total flux, using the central difference (CD) scheme.
     */
    void discretize_cd(Matrix& mat) const
    {
        discretize_cd(m_C,m_D,mat);
    }
    static constexpr double bernoulli(double x)
    {
        return std::abs(x) < 1e-9 ? 1-x/2 : x/std::expm1(x);
    }
    /** \a k is a face index
     *  NOTE: this requires that update has been called, so the Peclet number
     *  is up-to-date.
     */
    FluxCoefs calc_coefs_sg_contrib(const Vector& C, const Vector& D, Grid::Index k) const
    {
        const double du_Lf = m_grid.duCell(k-1) / 2;
        const double du_LH = m_grid.duNode(k);
        const double cf = du_Lf / du_LH;
        const double Pe = m_Pe(k);
        const double ber = bernoulli(Pe);
        /* NOTE: for small P we use the approximation
         * expm1(Pc)/expm1(P) = (1 + Pc + ... -1)/(1 + P + ... -1) = c + ...
         */
	const double expm1_ratio = std::abs(Pe)<1e-9 ? cf : std::expm1(Pe*cf)/std::expm1(Pe);
        const double A = D[k]/du_LH*ber*std::exp(Pe*cf) - C[k]*expm1_ratio;
        const double B = A + C[k];
        return { A, B };
    }
    /** Flux contribution, using the Scharfetter-Gummel (exponential) scheme.
     *  This uses the contribution's C and D, and the Peclet number that is
     *  based on the sums of the convective and diffusive contributions.
     */
    void discretize_sg(const Vector& C, const Vector& D, Matrix& mat) const
    {
        mat.fill(0.0);
        for (Grid::Index k = 1; k < m_grid.nCells(); ++k)
        {
            mat.coeffRef(k, k) = 0;

            const double du_we = m_grid.duCell(k);
            if (k > 0)
            {
                const auto coefs = calc_coefs_sg_contrib(C,D,k);
                mat.coeffRef(k, k - 1)  = +coefs.B / du_we;
                mat.coeffRef(k, k)     += -coefs.A / du_we;
            }

            if (k < m_grid.nCells() - 1)
            {
                const auto coefs = calc_coefs_sg_contrib(C,D,k+1);
                mat.coeffRef(k, k)     += -coefs.B / du_we;
                mat.coeffRef(k, k + 1)  = +coefs.A / du_we;
            }
        }
    }
    /** \a k is a face index
     *  NOTE: this requires that update has been called, so the members
     *  m_C, m_D and m_Pe are up-to-date.
     *  This has the same effect as calling the three-argument overload with
     *  members m_C and m_D, but this special case allows some simplifications.
     */
    FluxCoefs calc_coefs_sg(Grid::Index k) const
    {
        const double du_WP = m_grid.duNode(k);
        const double ber = bernoulli(m_Pe(k));
        const double A = m_D[k]/du_WP*ber;
        const double B = m_D[k]/du_WP*(ber+m_Pe(k));
        return { A, B };
    }
    /** Discretize the total flux, using the Scharfetter-Gummel (exponential)
     *  scheme.
     */
    void discretize_sg(Matrix& mat) const
    {
        mat.fill(0.0);
        for (Grid::Index k = 1; k < m_grid.nCells(); ++k)
        {
            mat.coeffRef(k, k) = 0;

            const double du_we = m_grid.duCell(k);
            if (k > 0)
            {
                const auto coefs = calc_coefs_sg(k);
                mat.coeffRef(k, k - 1)  = +coefs.B / du_we;
                mat.coeffRef(k, k)     += -coefs.A / du_we;
            }

            if (k < m_grid.nCells() - 1)
            {
                const auto coefs = calc_coefs_sg(k+1);
                mat.coeffRef(k, k)     += -coefs.B / du_we;
                mat.coeffRef(k, k + 1)  = +coefs.A / du_we;
            }
        }
    }
private:
    const Grid& m_grid;
    Contributions m_contribs;
    Vector m_C;
    Vector m_D;
    Vector m_Pe;
};

} // name loki

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

    std::ofstream osf("f_" + suffix + ".dat");
    for (loki::Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        osf << grid.getCells()[k] << '\t' << eedf[k] << '\t' << eedf_an[k] << std::endl;
    }
    return eedf;
}

int main()
{
    using namespace loki;
    //const double eps = std::numeric_limits<double>::epsilon();

    constexpr Grid::Index Nc = 1001;
    constexpr double u_max = 2e-19/Constant::electronCharge;
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

        solveCase(grid,mat,"old");
    }

    /* the following calculations solve the same problem, but using the
     * ConvectionDiffusionDiscretizer class. First we register the (same)
     * elastic and field operators, next we call updateCD on both to update
     * their convection and diffusion coefficients. Then we call update()
     * on the discretizer object, which will update the summed convection
     * and diffusion coefficients and the grid Peclet number.
     */

    ConvectionDiffusionDiscretizer disc(grid);
    disc.register_term(op_elast);
    disc.register_term(op_field);

    op_elast.updateCD(grid,sigma_f,Tg);
    op_field.updateCD(grid,sigma_f,EoN,WoN,CIEff);
    disc.update();

    Matrix mat(grid.nCells(),grid.nCells());

    /* 2. Use ConvectionDiffusionDiscretizer to do the central difference
     *    discretization of the total flux.
     */

    mat.fill(0.0);
    disc.discretize_cd(mat);
    solveCase(grid,mat,"cd");

    /* 3. Use ConvectionDiffusionDiscretizer to do the Scharfetter-Gummel
     *    discretization of the total flux.
     */

    mat.fill(0.0);
    disc.discretize_sg(mat);
    solveCase(grid,mat,"sg");

    /* 4. Use ConvectionDiffusionDiscretizer to do the Scharfetter-Gummel
     *    discretization of the individual flux contributions. Aggregate
     *    the results afterwards.
     */

    Matrix mat_elast(grid.nCells(),grid.nCells());
    mat_elast.fill(0.0);
    disc.discretize_sg(op_elast.conv_coeff(),op_elast.diff_coeff(),mat_elast);

    Matrix mat_field(grid.nCells(),grid.nCells());
    mat_field.fill(0.0);
    disc.discretize_sg(op_field.conv_coeff(),op_field.diff_coeff(),mat_field);

    mat = mat_elast + mat_field;
    const auto eedf = solveCase(grid,mat,"sg_sep");

    /* NOTE: at this point, the first rows of the matrices are modified to
     * incorporate the normalization constraint. But we do not care, since
     * the row is multiplied with a zero energy value in the calculation of
     * the power terms using the expressions below.
     *
     * This statement of the power terms is explained in the cpp_notes document.
     */
    const Vector gamma_u_du = SI::gamma*grid.duCells().array()*grid.getCells().array();
    const double P_elast = gamma_u_du.dot(mat_elast*eedf);
    const double P_field = gamma_u_du.dot(mat_field*eedf);

    std::cout << "Power balance --- volumetric power/(n_e*N) in eV*m^3/s:" << std::endl;
    std::cout << " Elastic: " << std::showpos << P_elast << std::endl;
    std::cout << " Field:   " << std::showpos << P_field << std::endl;
    std::cout << " Net:     " << std::showpos << (P_elast+P_field) << std::endl;
}
