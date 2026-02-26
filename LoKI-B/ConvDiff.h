/** \file
 *
 *  Interfaces of classes for managing and discretizing convective-diffusive
 *  fluxes and their divergence in energy space.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2025 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
 *  M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
 *  L. L. Alves
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  \author Daan Boer and Jan van Dijk
 *  \date   February 2026
 */

#ifndef LOKI_CPP_CONVDIFF_H
#define LOKI_CPP_CONVDIFF_H

#include "LoKI-B/Grid.h"
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

namespace loki {

/** Base class for convective-diffusive flux terms. For such terms, the flux
 *  (or flux contribution) is modelled as H = conv*eedf - diff*d eedf/du.
 *  Members conv_coeff and diff_coeff represent the face values of the
 *  convection and diffusion coefficients.
 *
 *  This class only provides read-only access to these coefficients. Derived
 *  classes that implement particular convective-diffusive terms can modify
 *  the (protected) members m_conv_coeff and m_diff_coeff.
 *
 *  \seealso ConvectionDiffusionTerms
 *
 *  \author Daan Boer and Jan van Dijk
 *  \date   February 2026
 */
class ConvectionDiffusionOperator
{
  public:
    const Vector &conv_coeff() const { return m_conv_coeff; }
    const Vector &diff_coeff() const { return m_diff_coeff; }
  protected:
    Vector m_conv_coeff;
    Vector m_diff_coeff;
};

/** This convective-diffusive term manages a std::vector of pointers to
 *  external convective-diffusive operators. It calculates the base class'
 *  convection and diffusive coefficients as the sums of the contributions.
 *  These sums are re-calculated when member update() is called. Terms can
 *  be added with member register_term.
 *
 *  \seealso ConvectionDiffusionOperator
 *
 *  \author Jan van Dijk
 *  \date February 2026
 */
class ConvectionDiffusionTerms : public ConvectionDiffusionOperator
{
public:
    /// Type of the container of pointers to the drift-diffusive terms
    using Terms = std::vector<const ConvectionDiffusionOperator*>;
    /// Type of the size of the Terms contains
    using size_type = Terms::size_type;
    /// returns a constant reference to the drift-diffusive terms
    const Terms& terms() const { return m_terms; }
    /** Register the pointer to \a term with the term container.
     *  The lifetime of \a term should exceed that of the present class.
     */
    void register_term(const ConvectionDiffusionOperator& term)
    {
        m_terms.push_back(&term);
    }
    /** Re-calculate the sum of the convection and diffusion coefficients of
     *  the terms and store the results in the inherited members m_conv_coeff
     *  and m_diff_coeff. Those are resized, if necessary. The coefficient
     *  vectors of the contributions must be all the same, otherwise the
     *  bahaviour is undefined (in debug mode, an assert will be triggered).
     *  If the term container is empty, members m_conv_coeff and m_diff_coeff
     *  will be resized to size 0.
     */
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
    /// The container of pointers to the drift-diffusive terms.
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
     *  Gamma_f = B*f_B - A*f_A, where f_B and f_A are the field values before
     *  and after the face. For a face with index n, the flux expression is
     *  Gamma[n] = B*f[n-1]+A*F[n] (see the Grid class for an explanation of the
     *  layout and numbering of faces and cells). If either field value does not
     *  exist (coefficient A or B at the upper and lower energy grid boundary,
     *  respectively), the corresponding coefficient is set to a (quiet) NaN by
     *  the functions in this namespace.
     */
    struct LocalFluxCoefficients
    {
        /// coefficient of the field value in the cell after the face
        double A;
        /// coefficient of the field value in the cell before the face
        double B;
    };
    /** This structure provides member calc_coefs for the central difference
     *  scheme.
     */
    struct CD
    {
        /** Calculate the flux coefficients A and B for the flux at face \a k
         *  using the central difference scheme.
         *  At the upper grid boundary, a ghost cell point is introduced
         *  at u=u_max, and the requirement that the *total* flux (calculated
         *  from \a sum, and using the central difference scheme as well) is
         *  zero is used to eliminate the ghost point contribution, using a
         *  modified coefficient B. (Coefficient A should not be used, is value
         *  is set to NaN.) At the lower boundary, coefficient B is set to Nan,
         *  and B=0.0: there the flux is zero.
         */
        static LocalFluxCoefficients calc_coefs(
            const Grid& grid,
            const ConvectionDiffusionOperator& term,
            const ConvectionDiffusionOperator& sum,
            Grid::Index k)
        {
            if (k==0)
            {
                return { 0.0, std::numeric_limits<double>::quiet_NaN() };
            }
            /* note: at u=u_max (the upper boundary face) this results in
             * du_LF==du_LH, which amounts to placing the ghost point (see
             * below) at location u=u_max.
             */
            const double du_Lf = grid.duCell(k-1) / 2;
            const double du_LH = grid.duNode(k);
            const double cf = du_Lf / du_LH;
            const double D = term.diff_coeff()[k];
            const double C = term.conv_coeff()[k];
            const double A = D/du_LH - C*cf;
            const double B = A + C;
            if (k == grid.getNodes().size() - 1)
            {
                /* At the upper boundary, the total flux (represented by 'sum')
                 * is given by Bsum*f_B - Asum*f_g = 0, where f_g is the value
                 * in the ghost point, which is located on the grid boundary.
                 * Then
                 * Gamma = B*f_B - A*f_g = B*f_B - A*f_B*(Bsum/Asum) = (B - A*(Bsum/Asum))*f_B.
                 */
                const double Dsum = sum.diff_coeff()[k];
                const double Csum = sum.conv_coeff()[k];
                const double Asum = Dsum/du_LH - Csum*cf;
                const double Bsum = Asum + Csum;
                return { std::numeric_limits<double>::quiet_NaN(), B - A*Bsum/Asum };
            }
            return { A, B };
        }
    };
    /** This structure provides member calc_coefs for the Scharfetter-Gummel
     *  ('exponential') scheme.
     */
    struct SG
    {
        /** Calculate and return the bernoulli function for argument value \a x.
         *  This is defined as B(x)=x/(exp(x)-1) for x not equal to 0, and
         *  special value B(0)=lim_{x->0} B(x)=1. For small argument values
         *  (|x|<1e-9), the approximation 1-x/2 of the Taylor expansion
         *  x/(1+x+x^2/2+... -1) is used, for larger values of |x|, std::expm1
         *  is used for the numerator exp(x)-1.
         */
        static constexpr double bernoulli(double x)
        {
            return std::abs(x) < 1e-9 ? 1-x/2 : x/std::expm1(x);
        }
        /** Calculate the flux coefficients A and B for the flux at face \a k
         *  using the Scharfetter-Gummel scheme.
         *  At the upper grid boundary, a ghost cell point is introduced
         *  at u=u_max, and the requirement that the *total* flux (calculated
         *  from \a sum, and using theScharfetter-Gummel scheme as well) is
         *  zero is used to eliminate the ghost point contribution, using a
         *  modified coefficient B. (Coefficient A should not be used, is value
         *  is set to NaN.) At the lower boundary, coefficient B is set to Nan,
         *  and B=0.0: there the flux is zero.
         *
         *  The function is optimized for the case that \a term and \a sum are
         *  identical (more precisely: that they refer to the same object). In
         *  that case, the usual exponential scheme can be used, and the
         *  boundary fluxes are zero.
         */
        static LocalFluxCoefficients calc_coefs(
            const Grid& grid,
            const ConvectionDiffusionOperator& term,
            const ConvectionDiffusionOperator& sum,
            Grid::Index k)
        {
            if (k==0)
            {
                return { 0.0, std::numeric_limits<double>::quiet_NaN() };
            }
            const double du_LH = grid.duNode(k);
            const double Dsum = sum.diff_coeff()[k];
            const double Csum = sum.conv_coeff()[k];
            const double Pe = Csum*du_LH/Dsum;
            const double ber = bernoulli(Pe);
            if (&term==&sum)
            {
                /* This means that we are discretizing the total flux. This
                 * allows some simplifications of the evaluation, but the result
                 * should be the same as when the general code is used (barring
                 * round-off errors).
                 */
                if (k == grid.getNodes().size() - 1)
                {
                    return { std::numeric_limits<double>::quiet_NaN(), 0.0 };
                }
                const double A = Dsum/du_LH*ber;
                const double B = A + Csum;
                return { A, B };
            }
            // Handle the general case: term is a partial flux.
            const double du_Lf = grid.duCell(k-1) / 2;
            const double cf = du_Lf / du_LH;
            /* NOTE: for small P we use the approximation
             * expm1(Pc)/expm1(P) = (1 + Pc + ... -1)/(1 + P + ... -1) = c + ...
             */
            const double expm1_ratio = std::abs(Pe)<1e-9 ? cf : std::expm1(Pe*cf)/std::expm1(Pe);
            const double D = term.diff_coeff()[k];
            const double C = term.conv_coeff()[k];
            const double A = D/du_LH*ber*std::exp(Pe*cf) - C*expm1_ratio;
            const double B = A + C;
            if (k == grid.getNodes().size() - 1)
            {
                /* At the upper boundary, the total flux (represented by 'sum')
                 * is given by Bsum*f_B - Asum*f_g = 0, where f_g is the value
                 * in the ghost point, which is located on the grid boundary.
                 * Then
                 * Gamma = B*f_B - A*f_g = B*f_B - A*f_B*(Bsum/Asum) = (B - A*(Bsum/Asum))*f_B.
                 */
                const double Asum = Dsum/du_LH*ber;
                const double Bsum = Asum + Csum;
                return { std::numeric_limits<double>::quiet_NaN(), B - A*Bsum/Asum };
            }
            return { A, B };
        }
    };

} // namespace Schemes

/** Discretize the divergence of flux contribution \a term, using the SchemeTraits.
 *  This uses the \a term's convection and diffusion coefficients,
 *  and the Peclet number that is based on the \a sum of the convective and
 *  diffusive contributions.
 *
 *  \author Daan Boer and Jan van Dijk
 *  \date   February 2026
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
         * upper boundary, Gamma = B*f_B instead (calc_coefs modifies
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

/** Discretize the divergence of the (total) flux \a sum, using the SchemeTraits.
 *
 *  \author Jan van Dijk
 *  \date   February 2026
 */
template <class SchemeTraits>
void discretize_dflux_du(Matrix& mat, const Grid& grid, const ConvectionDiffusionOperator& sum)
{
    discretize_dflux_du<SchemeTraits>(mat,grid,sum,sum);
}

/** Evaluate flux contribution \a term, flux, using the SchemeTraits.
 *
 *  \author Jan van Dijk
 *  \date   February 2026
 */
template <class SchemeTraits>
void evaluate_flux_density(Vector& flux, const Grid& grid, const Vector& eedf, const ConvectionDiffusionOperator& term, const ConvectionDiffusionOperator& sum)
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
    assert(std::isnan(coefs.A));
#else
    flux[k] = 0.0;
#endif
}

/** Evaluate (total) flux \a sum, flux, using the SchemeTraits.
 *
 *  \author Jan van Dijk
 *  \date   February 2026
 */
template <class SchemeTraits>
void evaluate_flux_density(Vector& flux, const Grid& grid, const Vector& eedf, const ConvectionDiffusionOperator& sum)
{
    evaluate_flux_density<SchemeTraits>(flux,grid,eedf,sum,sum);
}

} // namespace loki

#endif // LOKI_CPP_CONVDIFF_H
