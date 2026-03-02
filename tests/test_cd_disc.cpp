/** \file
 *
 *  Unit tests for the ConvectionDiffusionOperator
 *
 *  \author Jan van Dijk
 *  \date   February 2026
 */

#include <fstream>
#include <ios>

#include "LoKI-B/Operators.h"
#include "LoKI-B/EedfUtilities.h"
/// \todo This is not used, but is should be: transform this code into a regtest.
#include "tests/TestUtilities.h"

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
void testScheme(const loki::Grid& grid, const loki::ConvectionDiffusionTerms& conv_diff_terms, const std::string& scheme_name)
{
    using namespace loki;

    Matrix mat(grid.nCells(),grid.nCells());

    /* 1. Do the discretization of the total flux, solve the equation and
     * print the results.
     */

    mat.fill(0.0);
    discretize_dflux_du<SchemeT>(mat,grid,conv_diff_terms);
    solveCase(grid,mat,scheme_name);

    /* 2. Do the discretization of the individual flux contributions, store the
     * matrices. Aggregate the results afterwards to discretize the total flux,
     * solve the equation, then evaluate and print the individual flux contributions
     * and the power terms.
     */

    std::vector<Matrix> mat_contribs(conv_diff_terms.terms().size());
    Matrix mat_term(grid.nCells(),grid.nCells());
    for (typename ConvectionDiffusionTerms::size_type t=0; t!=conv_diff_terms.terms().size(); ++t)
    {
        mat_contribs[t].resize(grid.nCells(),grid.nCells());
        mat_contribs[t].fill(0.0);
        discretize_dflux_du<SchemeT>(mat_contribs[t],grid,*conv_diff_terms.terms()[t],conv_diff_terms);
        mat += mat_contribs[t];
    }
    const auto eedf = solveCase(grid,mat,scheme_name+"_sep");

    std::vector<Vector> flux_contribs(conv_diff_terms.terms().size());
    Vector flux_sum(grid.getNodes().size());
    flux_sum.fill(0.0);
    for (typename ConvectionDiffusionTerms::size_type t=0; t!=conv_diff_terms.terms().size(); ++t)
    {
        flux_contribs[t].resize(grid.getNodes().size());
        evaluate_flux_density<SchemeT>(flux_contribs[t],grid,eedf,*conv_diff_terms.terms()[t],conv_diff_terms);
        flux_sum += flux_contribs[t];
    }
    Vector flux_total(grid.getNodes().size());
    evaluate_flux_density<SchemeT>(flux_total,grid,eedf,conv_diff_terms);
    std::ofstream ofs("flux_"+scheme_name+"_sep.dat");
    for (loki::Grid::Index k = 0; k < grid.getNodes().size(); ++k)
    {
        ofs << grid.getNodes()[k];
        for (const auto& f : flux_contribs)
        {
            ofs << '\t' << f[k];
        }
        ofs << '\t' << flux_sum[k] << '\t' << flux_total[k] << std::endl;
    }

    /* NOTE: the first rows of the matrices are not set up by the discretizers
     * at present. We do not care, since this row is multiplied with a zero
     * energy value in the calculation of the power terms using the expressions
     * below.
     *
     * This statement of the power terms is explained in the cpp_notes document.
     */
    const Vector u_du = grid.duCells().array()*grid.getCells().array();
    std::cout << "Power balance --- volumetric power/(n_e*N) in eV*m^3/s:" << std::endl;
    double P_sum=0;
    for (typename ConvectionDiffusionTerms::size_type t=0; t!=conv_diff_terms.terms().size(); ++t)
    {
        const double P_term = u_du.dot(mat_contribs[t]*eedf);
        std::cout << "Term #" << t << ": " << std::showpos << P_term << std::endl;
        P_sum += P_term;
    }
    std::cout << "    Sum: " << std::showpos << P_sum << std::endl;

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
    conv_diff_terms.registerTerm(op_elast);
    conv_diff_terms.registerTerm(op_field);

    op_elast.updateCD(grid,sigma_f,Tg);
    op_field.updateCD(grid,sigma_f,EoN,WoN,CIEff);

    conv_diff_terms.update();

    testScheme<Schemes::CD>(grid, conv_diff_terms, "cd");
    testScheme<Schemes::SG>(grid, conv_diff_terms, "sg");

    return 0;
}
