/** \file
 *
 *  Implementations of classes for managing and discretizing convective-diffusive
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

#include "LoKI-B/ConvDiff.h"
#include <cassert>
#include <limits>

/** When 0, (partial) boundary fluxes are the upper grid boundary are
 *  assumed to be zero. Otherwise, those fluxes will be evaluated using
 *  a ghost point on the upper boundary under the assumption that the
 *  *total* convective-diffusive flux is 0.
 */
#define LOKI_DISCRETIZE_BOUNDARY_FLUX 1

namespace loki {

void ConvectionDiffusionTerms::register_term(const ConvectionDiffusionOperator& term)
{
    m_terms.push_back(&term);
}

void ConvectionDiffusionTerms::update()
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

namespace Schemes {

LocalFluxCoefficients CD::calc_coefs(
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
        if (&term==&sum)
        {
            return { std::numeric_limits<double>::quiet_NaN(), 0.0 };
        }
#if LOKI_DISCRETIZE_BOUNDARY_FLUX
        /* At the upper boundary, the total flux (represented by 'sum')
         * is given by Bsum*f_B - Asum*f_g = 0, where f_g is the value
         * in the ghost point, which is located on the grid boundary. Then
         * Gamma = B*f_B - A*f_g = B*f_B - A*f_B*(Bsum/Asum) = (B - A*(Bsum/Asum))*f_B.
         */
        const double Dsum = sum.diff_coeff()[k];
        const double Csum = sum.conv_coeff()[k];
        const double Asum = Dsum/du_LH - Csum*cf;
        const double Bsum = Asum + Csum;
        return { std::numeric_limits<double>::quiet_NaN(), B - A*Bsum/Asum };
#else
        return { std::numeric_limits<double>::quiet_NaN(), 0.0 };
#endif
    }
    return { A, B };
}

LocalFluxCoefficients SG::calc_coefs(
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
#if LOKI_DISCRETIZE_BOUNDARY_FLUX
        /* At the upper boundary, the total flux (represented by 'sum')
         * is given by Bsum*f_B - Asum*f_g = 0, where f_g is the value
         * in the ghost point, which is located on the grid boundary. Then
         * Gamma = B*f_B - A*f_g = B*f_B - A*f_B*(Bsum/Asum) = (B - A*(Bsum/Asum))*f_B.
         */
        const double Asum = Dsum/du_LH*ber;
        const double Bsum = Asum + Csum;
        return { std::numeric_limits<double>::quiet_NaN(), B - A*Bsum/Asum };
#else
        return { std::numeric_limits<double>::quiet_NaN(), 0.0 };
#endif
    }
    return { A, B };
}

} // namespace Schemes

template <class SchemeTraits>
void discretize_dflux_du(
    Matrix& mat,
    const Grid& grid,
    const ConvectionDiffusionOperator& term,
    const ConvectionDiffusionOperator& sum)
{
    mat.fill(0.0);
    /** \todo We visit every face twice. This can be fixed by changing
     *  this in a loop over the faces, and handle the contributions to the
     *  equations in both the lower and upper adjacent cells.
     */
    for (Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        const double du_we = grid.duCell(k);
        mat.coeffRef(k, k) = 0;

        // contribution from the flux at the lower face of k:
        {
            const auto coefs = SchemeTraits::calc_coefs(grid,term,sum,k);
            if (k > 0)
            {
                mat.coeffRef(k, k - 1) = +coefs.B / du_we;
            }
            mat.coeffRef(k, k) += -coefs.A / du_we;
        }

        // contribution from the flux at the upper face of k:
        {
            const auto coefs = SchemeTraits::calc_coefs(grid,term,sum,k+1);
            if (k < grid.nCells()-1)
            {
                mat.coeffRef(k, k + 1) = +coefs.A / du_we;
            }
            mat.coeffRef(k, k) += -coefs.B / du_we;
        }
    }
}

template <class SchemeTraits>
void discretize_dflux_du(
    Matrix& mat, 
    const Grid& grid, 
    const ConvectionDiffusionOperator& sum)
{
    discretize_dflux_du<SchemeTraits>(mat,grid,sum,sum);
}

template <class SchemeTraits>
void evaluate_flux_density(
    Vector& flux,
    const Grid& grid,
    const Vector& eedf,
    const ConvectionDiffusionOperator& term,
    const ConvectionDiffusionOperator& sum)
{
    // first cell:
    {
        const Grid::Index k = 0;
        const auto coefs = SchemeTraits::calc_coefs(grid,term,sum,k);
        flux[k] = coefs.A*eedf[k];
        assert(std::isnan(coefs.B));
    }

    // internal cells:
    for (Grid::Index k = 1; k < grid.getNodes().size()-1; ++k)
    {
            const auto coefs = SchemeTraits::calc_coefs(grid,term,sum,k);
            flux[k] = coefs.B*eedf[k-1] - coefs.A*eedf[k];
    }

    // last cell:
    {
        const Grid::Index k = grid.getNodes().size()-1;
        const auto coefs = SchemeTraits::calc_coefs(grid,term,sum,k);
        flux[k] = coefs.B*eedf[k-1];
        assert(std::isnan(coefs.A));
    }
}

template <class SchemeTraits>
void evaluate_flux_density(
    Vector& flux, 
    const Grid& grid, 
    const Vector& eedf, 
    const ConvectionDiffusionOperator& sum)
{
    evaluate_flux_density<SchemeTraits>(flux,grid,eedf,sum,sum);
}


// explicit instantiations for Schemes::CD:

template void discretize_dflux_du<Schemes::CD>(
    Matrix& mat, 
    const Grid& grid, 
    const ConvectionDiffusionOperator& term, 
    const ConvectionDiffusionOperator& sum);

template void discretize_dflux_du<Schemes::CD>(
    Matrix& mat, 
    const Grid& grid, 
    const ConvectionDiffusionOperator& sum);

template void evaluate_flux_density<Schemes::CD>(
    Vector& flux, 
    const Grid& grid, 
    const Vector& eedf, 
    const ConvectionDiffusionOperator& term, 
    const ConvectionDiffusionOperator& sum);

template void evaluate_flux_density<Schemes::CD>(
    Vector& flux, 
    const Grid& grid, 
    const Vector& eedf, 
    const ConvectionDiffusionOperator& sum);

// explicit instantiations for Schemes::SG:

template void discretize_dflux_du<Schemes::SG>(
    Matrix& mat, 
    const Grid& grid, 
    const ConvectionDiffusionOperator& term, 
    const ConvectionDiffusionOperator& sum);

template void discretize_dflux_du<Schemes::SG>(
    Matrix& mat, 
    const Grid& grid, 
    const ConvectionDiffusionOperator& sum);

template void evaluate_flux_density<Schemes::SG>(
    Vector& flux, 
    const Grid& grid, 
    const Vector& eedf, 
    const ConvectionDiffusionOperator& term, 
    const ConvectionDiffusionOperator& sum);

template void evaluate_flux_density<Schemes::SG>(
    Vector& flux, 
    const Grid& grid, 
    const Vector& eedf, 
    const ConvectionDiffusionOperator& sum);

} // namespace loki
