#include "LoKI-B/Grid.h"
#include "LoKI-B/EedfUtilities.h"
#include "LoKI-B/ElectronKinetics.h"
#include <vector>

namespace loki {

/** A convection-diffusion operator. Its members C and D are defined such that,
 *  on node i, the flux is given by H(u_i) = C_i*eedf(u_i) - D_i*(d eedf/du)_i
 *  (convection and diffusion coefficients).
 *
 *  The class provides members clear and operator+=(other) to facilitate
 *  forming a compound operator that contains summed C and D values.
 */
class ConvDiffOperator
{
public:
    ConvDiffOperator(const Grid& grid, const std::string& name)
    : m_grid(grid), m_name(name)
    {
        clear();
    }
    const Grid& grid() const { return m_grid; }
    const std::string& name() const { return m_name; }
    void clear()
    {
        m_C.resize(m_grid.getNodes().size());
        m_D.resize(m_grid.getNodes().size());
        m_C.fill(0.0);
        m_D.fill(0.0);
    }
    ConvDiffOperator& operator+=(const ConvDiffOperator& other)
    {
        m_C += other.m_C;
        m_D += other.m_D;
        return *this;
    }
private:
    const Grid& m_grid;
    const std::string m_name;
public:
    Vector m_C;
    Vector m_D;
};

/* The ConvDiffDiscretizer discretizes a convection-diffusion problem.
 * Multiple ConvDiffOperator contributions can be registered (member
 * register), the class also has a member m_cd_total that stores the
 * total (summed) C and D of the contributions and the grid Peclet
 * numbers at the faces, based on the total C and D.
 *
 * After registering contributions, or when either contribution's C
 * and/or D has changed, member update should be called to update the
 * total C and D and the Peclet numbers.
 *
 * Subsequently, member discretize can be called to *add* to the
 * Boltzmann matrix (if necessary, that should be cleared by the
 * user before calling this function).
 *
 * A second overload of discretize only discretizes one of the registered
 * contributions. Note that such contribution still depends on the global
 * Peclet number. It is possible to consider only the convective or
 * diffusive parts ('downflux/upflux') --- useful for power calculations.
 * Member calc_power calculates the power associated with such contribution
 * (very inefficient at present, it sets up a local matrix to carry out
 * the calculation).
 */
class ConvDiffDiscretizer
{
public:
    using ConvDiffOperators = std::vector<const ConvDiffOperator*>;
    ConvDiffDiscretizer(const Grid& grid)
    : m_cd_total(grid, "combined_flux")
    {
    }
    const Grid& grid() const { return m_cd_total.grid(); }
    void register_cd_operator(const ConvDiffOperator& op)
    {
        m_cd_ops.push_back(&op);
    }
    void update()
    {
        /* recalculate C and D of m_cd_total, representing the sums of the
         * contributions' C and D.
         */
        m_cd_total.clear();
        for (const auto& cd_op : m_cd_ops)
        {
            m_cd_total += *cd_op;
        }
        m_P = m_cd_total.m_C.cwiseProduct(grid().duNodes()).cwiseQuotient(m_cd_total.m_D);
    }
    /** For every cell [0,n), discretize the equation
     *  (H_{i+1/2}-H_{i_{i-1/2})/(du)_i. This is done by looping over the
     *  faces, rather than the cells to avoid visiting each face twice. The
     *  fluxes at the boundary faces are skipped (zero-flux boundary condition).
     */
    void discretize(SparseMatrix& mat) const
    {
        double A,B;
        // skip both boundary faces
        for (Grid::Index f=1; f!=grid().getNodes().size()-1; ++f)
        {
            calculateAB(A,B,f);
            // indices of cells before and after face f
            const Grid::Index fB = f-1;
            const Grid::Index fA = f;
            const double du = grid().duNode(f);
            if (fB!=0)
            {
                mat.coeffRef(fB,fB) += -B/du;
                mat.coeffRef(fB,fA) += +A/du;
            }
            if (fB!=grid().getNodes().size()-1)
            {
                mat.coeffRef(fA,fB) += +B/du;
                mat.coeffRef(fA,fA) += -A/du;
            }
        }
    }
    void discretize(SparseMatrix& mat, const ConvDiffOperator& cd_op, bool conv, bool diff) const
    {
        double Ax,Bx;
        // skip both boundary faces
        for (Grid::Index f=1; f!=grid().getNodes().size()-1; ++f)
        {
            calculateAB(Ax,Bx,f,cd_op,conv,diff);
            // indices of cells before and after face f
            const Grid::Index fB = f-1;
            const Grid::Index fA = f;
            const double du = grid().duNode(f);
            if (fB!=0)
            {
                mat.coeffRef(fB,fB) += -Bx/du;
                mat.coeffRef(fB,fA) += +Ax/du;
            }
            if (fB!=grid().getNodes().size()-1)
            {
                mat.coeffRef(fA,fB) += +Bx/du;
                mat.coeffRef(fA,fA) += -Ax/du;
            }
        }
#if 0
        /* use the ghost cell method to discretize the boundary flux for
         * this operator:
            H_bnd=0 <=> B*f_L - A*f_H = 0 <=> f_H = f_L*(B/A).
            Then
            H_x,bnd = B_x*f_L - A_x*f_H = B_x*f_L - A_x f_L*(B/A) = (B_x - A_x*(B/A))*f_L
         */
        /// \todo Check which du is (must be) used for calculating A,B, Ax, Bx.
        {
            Grid::Index f=grid().getNodes().size()-1;
            calculateAB(Ax,Bx,f,cd_op,conv,diff);
            double A,B;
            calculateAB(A,B,f);
            // index of cell before the boundary face
            const Grid::Index fB = f-1;
            const double du = grid().duNode(f);
            // std::cout << "MAT=" << mat.coeffRef(fB,fB) << ", DEL = " << (Bx-Ax*B/A)/du << std::endl;
            mat.coeffRef(fB,fB) += +(Bx-Ax*B/A)/du;
        }
#endif
    }
    const Vector& P() const { return m_P; }
    double calc_power(const ConvDiffOperator& op, const Vector& eedf, bool conv, bool diff) const
    {
        SparseMatrix mat(grid().nCells(),grid().nCells());
        discretize(mat,op,conv,diff);
        return grid().duCells().dot(grid().getCells().asDiagonal()*(mat*eedf));
    }
    const ConvDiffOperators& cd_ops() const { return m_cd_ops; }
    const ConvDiffOperator& cd_total() const { return m_cd_total; }
    void write_power(std::ostream& os, const ConvDiffOperator& term, const Vector& eedf) const
    {
        const double conv = calc_power(term, eedf, true, false);
        const double diff = calc_power(term, eedf, false, true);
        // note: this should be equal to conv + diff
        const double total = calc_power(term, eedf, true, true);
        os << "#  " << term.name() << std::endl;
        os << "#  " << std::right << std::setw(15) << "convective: " << conv << std::endl;
        os << "#  " << std::right << std::setw(15) << "diffusive:  " << std::showpos << diff << std::endl;
        os << "#  " << std::right << std::setw(15) << "total:      " << std::setw(30) << total << std::endl;
    }
private:
    double bernoulli(double P) const
    {
        /// \todo Use a threshold here, not 0?
        return P==0 ? 1.0 : P/std::expm1(P);
    }
    /** Calculate coefficients \a A and \a B such that the flux of a field
     *  at internal face \a f is given by H = B*field(f-1) - A*field(f).
     */
    void calculateAB(double& A, double& B, Grid::Index f) const
    {
        const double C=m_cd_total.m_C[f];
        const double D=m_cd_total.m_D[f];
        if (D==0)
        {
            A = C>=0 ? 0 : -C;
        }
        else
        {
            const double du = grid().duNode(f);
            const double P=m_P[f];
            // exponential scheme:
            const double Ber = bernoulli(P);
            A = (D/du)*Ber;
        }
        B = A + C;
    }
    /** Calculate coefficients \a A and \a B such that the flux of a field
     *  at internal face \a f due to convection-diffusion contribution \a cd_op
     *  is given by H = B*field(f-1) - A*field(f). The booleans \a conv and
     *  \a diff control if the convective and/or diffusive components are
     *  taken into account. This allows the discretization of only the
     *  'downflux' (convective) or 'upflux' (diffusive) parts of the flux,
     *  which is used in power balance calculations.
     */
    void calculateAB(double& A, double& B, Grid::Index f, const ConvDiffOperator& cd_op, bool conv, bool diff) const
    {
        A = 0.0;
        B = 0.0;
        // for the specific term:
        const double Cx=conv ? cd_op.m_C[f] : 0.0;
        const double Dx=diff ? cd_op.m_D[f] : 0.0;
        // for the combined terms:
        const double P=m_P[f];
        // grid data:
        const double du = grid().duNode(f);
        const double c = (grid().getNode(f)-grid().getCell(f-1))/du;
        if (Dx!=0)
        {
            /* if Dx!=0, also D!=0: P is finite
             * The expression is rewritten for positive P to avoid
             * an underflow-overflow combination for large positive P.
             */
            A += P>0 ? (Dx/du)*bernoulli(-P)*std::exp(P*(c-1)) : (Dx/du)*bernoulli(P)*std::exp(P*c);
        }
        if (Cx!=0)
        {
            using std::expm1;
            /* \todo Handle P=+/-inf: in principle, all Dx can be 0,
             * or verify that this is captured by the expression below.
             */
            A += c==0 ? -Cx*c : P>0 ? -Cx*std::exp(P*(c-1))*expm1(-P*c)/expm1(-P) : -Cx*expm1(P*c)/expm1(P);
        }
        B = A + Cx;
    }
    ConvDiffOperators m_cd_ops;
    ConvDiffOperator m_cd_total;
    Vector m_P;
};

} // namespace loki

int main()
{
    using namespace loki;
    Grid grid(1000,5);

    /* select one of the two testcases. true: maxwell at Te=Tg=0.3*eV
     * (no field). false: Druyvesteyn (Tg=0, E_N=5e-21*Vm^2).
     */
    const bool maxwell = false;

    const double mass_ratio = 1.0/1800;
    Vector sigma(grid.getNodes().size()); sigma.fill(1.0e-19);
    const double Tg_eV = maxwell ? 0.3 : 0;
    const double E_N = maxwell ? 0 : 5e-21;

    // set up an elastic operator
    ConvDiffOperator op_el(grid, "elastic");
    const Vector g_el = grid.getNodes().array()*grid.getNodes().array()*sigma.array()*2.0*mass_ratio;
    op_el.m_C = -g_el;
    op_el.m_D = g_el*Tg_eV;

    // set up a field operator
    ConvDiffOperator op_field(grid, "field");
    Vector g_E = E_N*E_N*(1./3)*grid.getNodes().array()/sigma.array();
    op_field.m_C.fill(0.0);
    op_field.m_D = g_E;

    /* set up a convection-diffusion (cd) discretizer. Register the elastic and
     * field operators.
     */
    ConvDiffDiscretizer cd(grid);
    cd.register_cd_operator(op_el);
    cd.register_cd_operator(op_field);

    /* By calling update, member m_cd_total of the discretizer is updated. Its
     * members m_C and m_D contain the sums of the members m_C and m_D of the
     * registered contributions. Subsequently, it updates the vector that holds
     * the grid Peclet numbers at the faces.
     */
    cd.update();

    /* Next, set up the 'Boltzmann matrix' M, defined such that dH/du is
     * approximated by M*eedf (eedf and 0 are column vectors).
     */
    SparseMatrix mat(grid.nCells(),grid.nCells());
    cd.discretize(mat);

    /* solveEEDF modifies its matrix argument. We let that operator on a copy
     * in case we want to use the unmodified 'mat' later on.
     */
    Matrix A(mat);
    Vector eedf(grid.nCells());
    loki::solveEEDF(eedf,A,grid);

    /* Calculate the mean energy and use that to construct the analytical
     * solution. Print lines containing f, f_anal and the relative error
     * to the console.
     */
    const double uMean = getMeanEnergy(eedf,grid);
    const Vector anal = makePrescribedEDF(grid,maxwell ? 1 : 2,uMean*2/3);
    for (Grid::Index i=0; i!=grid.nCells(); ++i)
    {
        std::cout << grid.getCell(i) << '\t' << eedf[i] << '\t' << anal[i] << '\t' << (eedf[i]-anal[i])/anal[i] << std::endl;
    }

    std::cout << "# power balance:" << std::endl;
    for (const auto& term : cd.cd_ops())
    {
        cd.write_power(std::cout,*term,eedf);
    }
    cd.write_power(std::cout,cd.cd_total(),eedf);

#if 0
    /* Some noisy additional tests. First we calculate the individual
     * contributions of the elastic and field operators to the boltzmann
     * matrix. These should add up to M.
     */

    SparseMatrix mat_el(grid.nCells(),grid.nCells());
    SparseMatrix mat_field(grid.nCells(),grid.nCells());
    cd.discretize(mat_el,op_el,true,true);
    cd.discretize(mat_field,op_field,true,true);

    // print M, M_el and M_field, as well as the error (M-M_el-M_field)
    std::cout << "M =\n" << mat << std::endl;
    std::cout << "M_el =\n" << mat_el << std::endl;
    std::cout << "M_field =\n" << mat_field << std::endl;
    std::cout << "Del = " << (mat-mat_el-mat_field) << std::endl;
#endif

}
