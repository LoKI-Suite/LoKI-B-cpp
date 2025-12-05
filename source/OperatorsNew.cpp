#include "LoKI-B/OperatorsNew.h"
#include "LoKI-B/Constant.h"
#include "LoKI-B/EedfCollisions.h"
#include "LoKI-B/EedfMixture.h"
#include "LoKI-B/Enumeration.h"
#include "LoKI-B/Grid.h"
#include "LoKI-B/GridOps.h"
#include "LoKI-B/Integrators.h"
#include "LoKI-B/Iterators.h"
#include "LoKI-B/Log.h"
#include <cmath>
#include <stdexcept>

namespace loki
{
namespace experimental
{
DriftDiffusionOperator::DriftDiffusionOperator(const Grid &grid)
    : drift_coeff(grid.getNodes().size()), diff_coeff(grid.getNodes().size())
{
}
const Vector &DriftDiffusionOperator::drift_coefficient()
{
    return this->drift_coeff;
}
const Vector &DriftDiffusionOperator::diffusion_coefficient()
{
    return this->diff_coeff;
}

ElasticOperator::ElasticOperator(const Grid &grid) : DriftDiffusionOperator(grid)
{
}
void ElasticOperator::evaluate(const Grid &grid, const Vector &elasticCrossSection, double T_gas)
{
    this->drift_coeff = 2. * grid.getNodes().cwiseAbs2().cwiseProduct(elasticCrossSection);
    this->diff_coeff = -(Constant::kBeV * T_gas) * this->drift_coefficient();
}

FieldOperator::FieldOperator(const Grid &grid) : DriftDiffusionOperator(grid)
{
}
void FieldOperator::evaluate(const Grid &grid, const Vector &total_cs, double EoN)
{
    this->drift_coeff.setZero();
    this->diff_coeff = -1. / 3. * (EoN * EoN) * grid.getNodes().cwiseQuotient(total_cs);
}

template <typename I>
void integrate_sink(const Grid &grid, const EedfCollision &col, const Eigen::VectorXd &eedf, Matrix &mat)
{
    const auto target_density = col.getTarget()->delta();

    const auto &cs = *col.crossSection;

    double u_start = 0.0;

    GridIterator grid_iter(grid, 0.0);
    InterpolatingIterator cs_int(cs.lookupTable().x(), cs.lookupTable().y(), u_start);
    InterpolatingIterator eedf_int(grid.getCells(), eedf, u_start);

    while (u_start < grid.uMax())
    {
        if (grid_iter.shouldAdvance(u_start))
        {
            grid_iter.advance();
        }
        if (cs_int.shouldAdvance(u_start))
        {
            cs_int.advance();
        }
        if (eedf_int.shouldAdvance(u_start))
        {
            eedf_int.advance();
        }

        double u_end = std::min({grid_iter.xHigh(), cs_int.switchOn(), eedf_int.switchOn()});

        I::setRows(u_start, u_end, 0.0, -target_density, grid_iter, cs_int, eedf_int, mat);

        u_start = u_end;
    }
}

template <typename I>
void integrate_source(const Grid &grid, const EedfCollision &col, const Vector &eedf, Matrix &mat)
{
    const auto target_density = col.getTarget()->delta();
    const auto &cs = *col.crossSection;

    double u_start = 0.0;
    const double u_offset = cs.threshold();

    GridIterator grid_iter(grid, 0.0);
    InterpolatingIterator cs_int(cs.lookupTable().x(), cs.lookupTable().y(), u_offset);
    InterpolatingIterator eedf_int(grid.getCells(), eedf, u_offset);

    while (u_start + u_offset < grid.uMax())
    {
        if (grid_iter.shouldAdvance(u_start))
        {
            grid_iter.advance();
        }
        if (cs_int.shouldAdvance(u_start))
        {
            cs_int.advance();
        }
        if (eedf_int.shouldAdvance(u_start))
        {
            eedf_int.advance();
        }

        double u_end = std::min({grid_iter.xHigh(), cs_int.switchOn(), eedf_int.switchOn()});

        I::setRows(u_start, u_end, u_offset, target_density, grid_iter, cs_int, eedf_int, mat);

        u_start = u_end;
    }
}

template <typename I>
void integrate_sup_sink(const Grid &grid, const EedfCollision &col, const Eigen::VectorXd &eedf, Matrix &mat)
{
    const double product_density = col.m_rhsHeavyStates[0]->delta();
    const double swRatio = col.getTarget()->statisticalWeight / col.m_rhsHeavyStates[0]->statisticalWeight;

    const auto &cs = *col.crossSection;

    double u_start = 0.0;
    const double u_offset = cs.threshold();

    GridIterator grid_iter(grid, 0.0);
    InterpolatingIterator cs_int(cs.lookupTable().x(), cs.lookupTable().y(), u_offset);
    InterpolatingIterator eedf_int(grid.getCells(), eedf, 0.0);

    while (u_start + u_offset < grid.uMax())
    {
        if (grid_iter.shouldAdvance(u_start))
        {
            grid_iter.advance();
        }
        if (cs_int.shouldAdvance(u_start))
        {
            cs_int.advance();
        }
        if (eedf_int.shouldAdvance(u_start))
        {
            eedf_int.advance();
        }

        double u_end = std::min({grid_iter.xHigh(), cs_int.switchOn(), eedf_int.switchOn()});

        I::setRows(u_start, u_end, u_offset, -product_density * swRatio, grid_iter, cs_int, eedf_int, mat);

        u_start = u_end;
    }
}

template <typename I>
void integrate_sup_source(const Grid &grid, const EedfCollision &col, const Vector &eedf, Matrix &mat)
{
    const double product_density = col.m_rhsHeavyStates[0]->delta();
    const double swRatio = col.getTarget()->statisticalWeight / col.m_rhsHeavyStates[0]->statisticalWeight;

    const auto &cs = *col.crossSection;

    double u_start = 0.0;
    const double u_offset = cs.threshold();

    GridIterator grid_iter(grid, u_offset);
    InterpolatingIterator cs_int(cs.lookupTable().x(), cs.lookupTable().y(), u_offset);
    InterpolatingIterator eedf_int(grid.getCells(), eedf, 0.0);

    while (u_start + u_offset < grid.uMax())
    {
        if (grid_iter.shouldAdvance(u_start + u_offset))
        {
            grid_iter.advance();
        }
        if (cs_int.shouldAdvance(u_start))
        {
            cs_int.advance();
        }
        if (eedf_int.shouldAdvance(u_start))
        {
            eedf_int.advance();
        }

        double u_end = std::min({grid_iter.xHigh() - u_offset, cs_int.switchOn(), eedf_int.switchOn()});

        if (u_start == u_end)
        {
            exit(1);
        }

        I::setRows(u_start, u_end, u_offset, product_density * swRatio, grid_iter, cs_int, eedf_int, mat);

        u_start = u_end;
    }
}

InelasticOperator::InelasticOperator(const Grid &grid)
    : inelasticMatrix(grid.nCells(), grid.nCells()), superelasticMatrix(grid.nCells(), grid.nCells())
{
}
void InelasticOperator::evaluate(const Grid &grid, const Vector &eedf, const EedfMixture &mixture)
{
    inelasticMatrix.setZero();

    for (const auto &cd : mixture.collision_data().data_per_gas())
    {
        for (auto vecIndex : {CollisionType::excitation, CollisionType::vibrational, CollisionType::rotational})
        {
            for (const auto &collision : cd.collisions(vecIndex))
            {
                if (collision->crossSection->threshold() >= grid.uMax())
                    continue;

                integrate_sink<LinIntegrator>(grid, *collision, eedf, this->inelasticMatrix);
                integrate_source<LinIntegrator>(grid, *collision, eedf, this->inelasticMatrix);

                if (collision->isReverse())
                {
                    integrate_sup_sink<LinIntegrator>(grid, *collision, eedf, this->inelasticMatrix);
                    integrate_sup_source<LinIntegrator>(grid, *collision, eedf, this->inelasticMatrix);
                }
            }
        }
    }

    // NOTE: This line is only required when the drift diffusion terms are also
    // divided by the cell width.
    // inelasticMatrix.array().rowwise() /= grid.duCells().array().transpose();
}

// TODO: Write equal sharing source integral using the new iterator/integrator classes.
// void equal_sharing_source(const Grid &grid, const Vector &eedf, Matrix &matrix, const EedfCollision &col)
// {
//     const double threshold = col.crossSection->threshold();
//     const double target_density = col.getTarget()->delta();

//     double u_start = grid.getNode(0) + threshold;

//     // The grid iterator iterates over the source cells of the grid.
//     Grid::Index i_grid = 0;

//     // Find the target cell in which the left face of the first source cell lands.
//     Grid::Index j_grid =
//         std::upper_bound(grid.getNodes().begin(), grid.getNodes().end(), u_start) - grid.getNodes().begin() - 1;

//     // Compute initial f_slope.
//     double f_slope = u_start >= grid.getCell(j_grid) || j_grid == 0
//                          ? std::log(eedf[j_grid + 1] / eedf[j_grid]) / grid.duNode(j_grid + 1)
//                          : std::log(eedf[j_grid] / eedf[j_grid - 1]) / grid.duNode(j_grid);
//     double u_f_start = u_start >= grid.getCell(j_grid) || j_grid == 0 ? grid.getCell(j_grid) : grid.getCell(j_grid -
//     1);

//     // Get the raw cross section data.
//     const auto &cs = col.crossSection->lookupTable().y();
//     const auto &cs_energy = col.crossSection->lookupTable().x();

//     Grid::Index j_cs = std::upper_bound(cs_energy.begin(), cs_energy.end(), u_start) - cs_energy.begin() - 1;

//     double sig_slope = (cs[j_cs + 1] - cs[j_cs]) / (cs_energy[j_cs + 1] - cs_energy[j_cs]);
//     double u_sig_start = cs_energy[j_cs];
//     double u_sig_next = cs_energy[j_cs + 1];
//     double sig_start = cs[j_cs];

//     while (u_start < grid.uMax())
//     {
//         // If the current integration domain is outside the current cross
//         // section cell, move to the next cell.
//         if (u_start >= u_sig_next)
//         {
//             j_cs++;

//             if (j_cs == cs.size() - 1)
//             {
//                 sig_slope = 0.;
//                 u_sig_start = 0.;
//                 u_sig_next = std::numeric_limits<double>::max();
//                 sig_start = 0.;
//             }
//             else
//             {
//                 sig_slope = (cs[j_cs + 1] - cs[j_cs]) / (cs_energy[j_cs + 1] - cs_energy[j_cs]);
//                 u_sig_start = cs_energy[j_cs];
//                 u_sig_next = cs_energy[j_cs + 1];
//                 sig_start = cs[j_cs];
//             }
//         }
//         // If the current integration domain is outside the current source cell,
//         // move to the next grid cell.
//         if (u_start >= 2.0 * grid.getNode(i_grid + 1) + threshold)
//         {
//             i_grid++;
//         }
//         // If the current integration domain is outside the current target cell,
//         // move to the next grid cell.
//         if (u_start >= grid.getNode(j_grid + 1))
//         {
//             j_grid++;
//             u_f_start = grid.getCell(j_grid);
//         }
//         // If the current integration domain starts at the current target cell
//         // center, recompute the slope of the eedf.
//         if (u_start == grid.getCell(j_grid) && j_grid < grid.nCells() - 1)
//         {
//             f_slope = std::log(eedf[j_grid + 1] / eedf[j_grid]) / grid.duNode(j_grid + 1);
//         }
//         // The end of the current integration domain is either the next source
//         // cell center, the next cross section entry, the next target cell face,
//         // or the grid boundary.
//         const double u_end =
//             std::min(std::min(u_sig_next, 2.0 * grid.getNode(i_grid + 1) + threshold),
//                      u_start >= grid.getCell(j_grid) ? grid.getNode(j_grid + 1) : grid.getCell(j_grid));

//         double integral_start = collision_integral(u_sig_start, sig_start, sig_slope, u_f_start, f_slope, u_start);
//         double integral_end = collision_integral(u_sig_start, sig_start, sig_slope, u_f_start, f_slope, u_end);

//         matrix(i_grid, j_grid) += 4.0 * target_density * (integral_end - integral_start);

//         u_start = u_end;
//     }
// }

IonizationOperator::IonizationOperator(const Grid &grid, IonizationOperatorType type)
    : operatorType(type), ionizationMatrix(grid.nCells(), grid.nCells())
{
}

// NOTE: For now this only implements conservative ionization and equal sharing.
void IonizationOperator::evaluate(const Grid &grid, const Vector &eedf, const EedfMixture &mixture)
{
    ionizationMatrix.setZero();

    for (const auto &cd : mixture.collision_data().data_per_gas())
    {
        for (const auto &collision : cd.collisions(CollisionType::ionization))
        {
            if (collision->crossSection->threshold() > grid.uMax())
                continue;

            integrate_sink<LinIntegrator>(grid, *collision, eedf, ionizationMatrix);

            switch (operatorType)
            {
            case IonizationOperatorType::conservative:
                integrate_source<LinIntegrator>(grid, *collision, eedf, ionizationMatrix);
                break;
            default:
                throw std::runtime_error("For now only conservative ionization is supported.");
                // case IonizationOperatorType::equalSharing:
                //     equal_sharing_source(grid, eedf, ionizationMatrix, *collision);
                //     break;
            }
        }
    }
}

SpatialGrowthOperator::SpatialGrowthOperator(const Grid &grid) : DriftDiffusionOperator(grid)
{
}

void SpatialGrowthOperator::evaluate(const Grid &grid, const Vector &eedf, const Vector &total_cs, double EoN,
                                     const Matrix &ionizationMatrix, const Matrix &attachmentMatrix)
{
    const auto coefsCI = SI::gamma * grid.duCells().transpose() * (ionizationMatrix + attachmentMatrix);
    const double nu_eff = eedf.dot(coefsCI);

    const auto total_cs_cells = (total_cs.head(total_cs.size() - 1) + total_cs.tail(total_cs.size() - 1)) / 2.;

    const Vector D0 = grid.getCells().array() / (3. * total_cs_cells).array();
    const Vector D0Nodes = grid.getNodes().array() / (3. * total_cs).array();

    // This is 33a from \cite Manual_1_0_0
    const double ND = SI::gamma * energyIntegral(grid, D0, eedf);
    const double muE = -SI::gamma * fNodegPrimeEnergyIntegral(grid, D0Nodes, eedf) * EoN;

    const double discriminant = muE * muE - 4 * nu_eff * ND;
    const double alphaRedEff = (discriminant < 0.) ? nu_eff / muE : (muE - std::sqrt(discriminant)) / (2 * ND);

    const auto shared_term = D0Nodes.array() * EoN * alphaRedEff;

    this->drift_coeff = -D0Nodes.array() * alphaRedEff * alphaRedEff - shared_term;
    this->diff_coeff = -shared_term;

    Log<Message>::Warning("Mobility: ", muE / EoN);
}
} // namespace experimental
} // namespace loki
