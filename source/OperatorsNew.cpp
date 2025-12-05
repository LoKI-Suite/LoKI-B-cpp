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

// Assumes linear interpolation for the cross sections, and logarithmic
// interpolation for the eedf.
double collision_integral(double u_sig_start, double sig_start, double sig_slope, double u_f_start, double f_slope,
                          double u)
{
    const double expr =
        (2 * sig_slope / std::pow(f_slope, 3.) + (u * (sig_start + (u - u_sig_start) * sig_slope)) / f_slope -
         (sig_start + (2 * u - u_sig_start) * sig_slope) / std::pow(f_slope, 2.));
    return std::exp(f_slope * (u - u_f_start)) * expr;
    // return std::expm1(f_slope * (u - u_f_start)) * expr + expr;
}

// Assumes linear interpolation for the cross sections, and logarithmic
// interpolation for the eedf. This function avoids a call to std::exp.
double collision_integral_adv(double u_sig_start, double sig_start, double sig_slope, double u_f_start, double f_div,
                              double delta_u_f, double u)
{
    // exp(ln(f1/f2) / (u2 - u1)) = exp(ln((f2/f1)^((u2 - u1)^-1)) = (f2/f1)^((u2-u1)^-1)
    const double f_slope = std::log(f_div) / delta_u_f;
    const double expr =
        (2 * sig_slope / std::pow(f_slope, 3.) + (u * (sig_start + (u - u_sig_start) * sig_slope)) / f_slope -
         (sig_start + (2 * u - u_sig_start) * sig_slope) / std::pow(f_slope, 2.));
    return std::pow(f_div, (u - u_f_start) / delta_u_f) * expr;
}

// Assumes linear interpolation for the cross sections, and a constant value for
// the eedf.
double collision_integral_simple(double u_sig_start, double sig_start, double sig_slope, double u_f_start,
                                 double f_slope, double u)
{
    return 0.5 * std::pow(u, 2.) * (sig_start - u_sig_start * sig_slope) + 1. / 3. * std::pow(u, 3.) * sig_slope;
}

// This version of the `collision_integral_sink` function directly interpolates
// cross sections from their raw data.
void collision_integral_sink(const Grid &grid, const Vector &eedf, Matrix &inelasticMatrix, const EedfCollision &col)
{
    const double target_density = col.getTarget()->delta();

    // The grid iterator iterates over the cells of the grid.
    Grid::Index i_grid = 0;

    // Compute initial f_slope using forward differences.
    double f_slope = std::log(eedf[1] / eedf[0]) / grid.duNode(1);
    double u_f_start = grid.getCell(0);

    // Get the raw cross section data.
    const auto &cs = col.crossSection->lookupTable().y();
    const auto &cs_energy = col.crossSection->lookupTable().x();

    // This initialization of the i_cs index allows to treat cross section that
    // do not start at 0 energy.
    Grid::Index i_cs = cs_energy[0] == 0. ? 0 : -1;

    double sig_slope = 0.;
    double u_sig_start = 0.;
    double u_sig_next = cs_energy[i_cs + 1];
    double sig_start = 0.;

    double u_start = 0;

    while (u_start < grid.uMax())
    {
        // If the current integration domain is outside the current cross
        // section cell, move to the next cell.
        if (u_start >= u_sig_next)
        {
            i_cs++;

            if (i_cs == cs.size() - 1)
            {
                sig_slope = 0.;
                u_sig_start = 0.;
                u_sig_next = std::numeric_limits<double>::max();
                sig_start = 0.;
            }
            else
            {
                sig_slope = (cs[i_cs + 1] - cs[i_cs]) / (cs_energy[i_cs + 1] - cs_energy[i_cs]);
                u_sig_start = cs_energy[i_cs];
                u_sig_next = cs_energy[i_cs + 1];
                sig_start = cs[i_cs];
            }
        }
        // If the current integration domain is outside the current grid cell,
        // move to the next grid cell.
        if (u_start >= grid.getNode(i_grid + 1))
        {
            i_grid++;
            u_f_start = grid.getCell(i_grid);
        }
        // If the current integration domain starts at the current cell center,
        // recompute the slope of the eedf.
        if (u_start == grid.getCell(i_grid) && i_grid < grid.nCells() - 1)
        {
            f_slope = std::log(eedf[i_grid + 1] / eedf[i_grid]) / grid.duNode(i_grid + 1);
        }

        // The end of the current integration domain is either the next cell
        // center, the next cross section entry, or the grid boundary.
        const double u_end = std::min(
            std::min(u_start >= grid.getCell(i_grid) ? grid.getNode(i_grid + 1) : grid.getCell(i_grid), u_sig_next),
            grid.uMax());

        double integral_start = collision_integral(u_sig_start, sig_start, sig_slope, u_f_start, f_slope, u_start);
        double integral_end = collision_integral(u_sig_start, sig_start, sig_slope, u_f_start, f_slope, u_end);

        inelasticMatrix(i_grid, i_grid) -= target_density * (integral_end - integral_start);

        u_start = u_end;
    }
}

// This version of the `collision_integral_source` function directly
// interpolates cross sections from their raw data.
void collision_integral_source(const Grid &grid, const Vector &eedf, Matrix &inelasticMatrix, const EedfCollision &col)
{
    const double threshold = col.crossSection->threshold();
    const double target_density = col.getTarget()->delta();

    double u_start = grid.getNode(0) + threshold;

    // The grid iterator iterates over the source cells of the grid.
    Grid::Index i_grid = 0;

    // Find the target cell in which the left face of the first source cell lands.
    Grid::Index j_grid =
        std::upper_bound(grid.getNodes().begin(), grid.getNodes().end(), u_start) - grid.getNodes().begin() - 1;

    // Compute initial f_slope.
    double f_slope = u_start >= grid.getCell(j_grid) || j_grid == 0
                         ? std::log(eedf[j_grid + 1] / eedf[j_grid]) / grid.duNode(j_grid + 1)
                         : std::log(eedf[j_grid] / eedf[j_grid - 1]) / grid.duNode(j_grid);
    double u_f_start = u_start >= grid.getCell(j_grid) || j_grid == 0 ? grid.getCell(j_grid) : grid.getCell(j_grid - 1);

    // Get the raw cross section data.
    const auto &cs = col.crossSection->lookupTable().y();
    const auto &cs_energy = col.crossSection->lookupTable().x();

    Grid::Index j_cs = std::upper_bound(cs_energy.begin(), cs_energy.end(), u_start) - cs_energy.begin() - 1;

    double sig_slope = (cs[j_cs + 1] - cs[j_cs]) / (cs_energy[j_cs + 1] - cs_energy[j_cs]);
    double u_sig_start = cs_energy[j_cs];
    double u_sig_next = cs_energy[j_cs + 1];
    double sig_start = cs[j_cs];

    while (u_start < grid.uMax())
    {
        // If the current integration domain is outside the current cross
        // section cell, move to the next cell.
        if (u_start >= u_sig_next)
        {
            j_cs++;

            if (j_cs == cs.size() - 1)
            {
                sig_slope = 0.;
                u_sig_start = 0.;
                u_sig_next = std::numeric_limits<double>::max();
                sig_start = 0.;
            }
            else
            {
                sig_slope = (cs[j_cs + 1] - cs[j_cs]) / (cs_energy[j_cs + 1] - cs_energy[j_cs]);
                u_sig_start = cs_energy[j_cs];
                u_sig_next = cs_energy[j_cs + 1];
                sig_start = cs[j_cs];
            }
        }
        // If the current integration domain is outside the current source cell,
        // move to the next grid cell.
        if (u_start >= grid.getNode(i_grid + 1) + threshold)
        {
            i_grid++;
        }
        // If the current integration domain is outside the current target cell,
        // move to the next grid cell.
        if (u_start >= grid.getNode(j_grid + 1))
        {
            j_grid++;
            u_f_start = grid.getCell(j_grid);
        }
        // If the current integration domain starts at the current target cell
        // center, recompute the slope of the eedf.
        if (u_start == grid.getCell(j_grid) && j_grid < grid.nCells() - 1)
        {
            f_slope = std::log(eedf[j_grid + 1] / eedf[j_grid]) / grid.duNode(j_grid + 1);
        }
        // The end of the current integration domain is either the next source
        // cell center, the next cross section entry, the next target cell face,
        // or the grid boundary.
        const double u_end =
            std::min(std::min(u_sig_next, grid.getNode(i_grid + 1) + threshold),
                     u_start >= grid.getCell(j_grid) ? grid.getNode(j_grid + 1) : grid.getCell(j_grid));

        double integral_start = collision_integral(u_sig_start, sig_start, sig_slope, u_f_start, f_slope, u_start);
        double integral_end = collision_integral(u_sig_start, sig_start, sig_slope, u_f_start, f_slope, u_end);

        inelasticMatrix(i_grid, j_grid) += target_density * (integral_end - integral_start);

        u_start = u_end;
    }
}

void collision_integral_sup_sink(const Grid &grid, const Vector &eedf, Matrix &inelasticMatrix,
                                 const EedfCollision &col)
{
    const double threshold = col.crossSection->threshold();

    const double product_density = col.m_rhsHeavyStates[0]->delta();
    const double swRatio = col.getTarget()->statisticalWeight / col.m_rhsHeavyStates[0]->statisticalWeight;

    double u_start = threshold;

    // The grid iterator iterates over the source cells of the grid.
    Grid::Index i_grid = 0;

    // Compute initial f_slope using forward differences.
    double f_slope = std::log(eedf[1] / eedf[0]) / grid.duNode(1);
    double u_f_start = grid.getCell(0) + threshold;

    // Find the target cell in which the left face of the first source cell lands.
    Grid::Index j_grid =
        std::upper_bound(grid.getNodes().begin(), grid.getNodes().end(), threshold) - grid.getNodes().begin() - 1;

    // Get the raw cross section data.
    const auto &cs = col.crossSection->lookupTable().y();
    const auto &cs_energy = col.crossSection->lookupTable().x();

    Grid::Index j_cs = std::upper_bound(cs_energy.begin(), cs_energy.end(), threshold) - cs_energy.begin() - 1;

    double sig_slope = (cs[j_cs + 1] - cs[j_cs]) / (cs_energy[j_cs + 1] - cs_energy[j_cs]);
    double u_sig_start = cs_energy[j_cs];
    double u_sig_next = cs_energy[j_cs + 1];
    double sig_start = cs[j_cs];

    while (u_start < grid.uMax())
    {
        // If the current integration domain is outside the current cross
        // section cell, move to the next cell.
        if (j_cs < cs.size() - 1 && u_start >= u_sig_next)
        {
            j_cs++;

            if (j_cs == cs.size() - 1)
            {
                sig_slope = 0.;
                u_sig_start = 0.;
                u_sig_next = std::numeric_limits<double>::max();
                sig_start = 0.;
            }
            else
            {
                sig_slope = (cs[j_cs + 1] - cs[j_cs]) / (cs_energy[j_cs + 1] - cs_energy[j_cs]);
                u_sig_start = cs_energy[j_cs];
                u_sig_next = cs_energy[j_cs + 1];
                sig_start = cs[j_cs];
            }
        }
        // If the current integration domain is outside the current source cell,
        // move to the next grid cell.
        if (u_start >= (grid.getNode(i_grid + 1) + threshold) * (1.0 - 1e-14))
        {
            i_grid++;
            u_f_start = grid.getCell(i_grid) + threshold;
        }
        // If the current integration domain is outside the current target cell,
        // move to the next grid cell.
        if (u_start >= grid.getNode(j_grid + 1))
        {
            j_grid++;
        }
        // If the current integration domain starts at the current target cell
        // center, recompute the slope of the eedf.
        if (u_start == grid.getCell(i_grid) + threshold && i_grid < grid.nCells() - 1)
        {
            f_slope = std::log(eedf[i_grid + 1] / eedf[i_grid]) / grid.duNode(i_grid + 1);
        }
        // The end of the current integration domain is either the next source
        // cell center, the next cross section entry, the next target cell face,
        // or the grid boundary.
        const double u_end =
            std::min(std::min(u_sig_next, grid.getNode(j_grid + 1)),
                     (u_start >= (grid.getCell(i_grid) + threshold) * (1.0 - 1e-14) ? grid.getNode(i_grid + 1)
                                                                                    : grid.getCell(i_grid)) +
                         threshold);

        double integral_start = collision_integral(u_sig_start, sig_start, sig_slope, u_f_start, f_slope, u_start);
        double integral_end = collision_integral(u_sig_start, sig_start, sig_slope, u_f_start, f_slope, u_end);

        inelasticMatrix(i_grid, i_grid) -= product_density * swRatio * (integral_end - integral_start);

        u_start = u_end;
    }
}

void collision_integral_sup_source(const Grid &grid, const Vector &eedf, Matrix &inelasticMatrix,
                                   const EedfCollision &col)
{
    const double threshold = col.crossSection->threshold();

    const double product_density = col.m_rhsHeavyStates[0]->delta();
    const double swRatio = col.getTarget()->statisticalWeight / col.m_rhsHeavyStates[0]->statisticalWeight;

    double u_start = threshold;

    // The grid iterator iterates over the source cells of the grid.
    Grid::Index i_grid = 0;

    // Compute initial f_slope using forward differences.
    double f_slope = std::log(eedf[1] / eedf[0]) / grid.duNode(1);
    double u_f_start = grid.getCell(0);

    // Find the target cell in which the left face of the first source cell lands.
    Grid::Index j_grid =
        std::upper_bound(grid.getNodes().begin(), grid.getNodes().end(), threshold) - grid.getNodes().begin() - 1;

    // Get the raw cross section data.
    const auto &cs = col.crossSection->lookupTable().y();
    const auto &cs_energy = col.crossSection->lookupTable().x();

    Grid::Index j_cs = std::upper_bound(cs_energy.begin(), cs_energy.end(), threshold) - cs_energy.begin() - 1;

    double sig_slope = (cs[j_cs + 1] - cs[j_cs]) / (cs_energy[j_cs + 1] - cs_energy[j_cs]);
    double u_sig_start = cs_energy[j_cs];
    double u_sig_next = cs_energy[j_cs + 1];
    double sig_start = cs[j_cs];

    while (u_start < grid.uMax())
    {
        // If the current integration domain is outside the current cross
        // section cell, move to the next cell.
        if (j_cs < cs.size() - 1 && u_start >= u_sig_next)
        {
            j_cs++;

            if (j_cs == cs.size() - 1)
            {
                sig_slope = 0.;
                u_sig_start = 0.;
                u_sig_next = std::numeric_limits<double>::max();
                sig_start = 0.;
            }
            else
            {
                sig_slope = (cs[j_cs + 1] - cs[j_cs]) / (cs_energy[j_cs + 1] - cs_energy[j_cs]);
                u_sig_start = cs_energy[j_cs];
                u_sig_next = cs_energy[j_cs + 1];
                sig_start = cs[j_cs];
            }
        }
        // If the current integration domain is outside the current source cell,
        // move to the next grid cell.
        if (u_start >= grid.getNode(i_grid + 1) + threshold)
        {
            i_grid++;
            u_f_start = grid.getCell(i_grid) + threshold;
        }
        // If the current integration domain is outside the current target cell,
        // move to the next grid cell.
        if (u_start >= grid.getNode(j_grid + 1))
        {
            j_grid++;
        }
        // If the current integration domain starts at the current target cell
        // center, recompute the slope of the eedf.
        if (u_start == grid.getCell(i_grid) + threshold && i_grid < grid.nCells() - 1)
        {
            f_slope = std::log(eedf[i_grid + 1] / eedf[i_grid]) / grid.duNode(i_grid + 1);
        }
        // The end of the current integration domain is either the next source
        // cell center, the next cross section entry, the next target cell face,
        // or the grid boundary.
        const double u_end =
            std::min(std::min(u_sig_next, grid.getNode(j_grid + 1)),
                     (u_start >= grid.getCell(i_grid) + threshold ? grid.getNode(i_grid + 1) : grid.getCell(i_grid)) +
                         threshold);

        double integral_start = collision_integral(u_sig_start, sig_start, sig_slope, u_f_start, f_slope, u_start);
        double integral_end = collision_integral(u_sig_start, sig_start, sig_slope, u_f_start, f_slope, u_end);

        inelasticMatrix(j_grid, i_grid) += product_density * swRatio * (integral_end - integral_start);

        u_start = u_end;
    }
}

template <typename I>
void integrate_sink(const Grid &grid, const EedfCollision &col, const Vector &eedf, Matrix &mat);

template <typename I>
void integrate_source(const Grid &grid, const EedfCollision &col, const Vector &eedf, Matrix &mat);

template <typename I>
void integrate_sup_sink(const Grid &grid, const EedfCollision &col, const Vector &eedf, Matrix &mat);

template <typename I>
void integrate_sup_source(const Grid &grid, const EedfCollision &col, const Vector &eedf, Matrix &mat);

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

                // collision_integral_sink(grid, eedf, this->inelasticMatrix, *collision);
                // collision_integral_source(grid, eedf, this->inelasticMatrix, *collision);
                integrate_sink<LinIntegrator>(grid, *collision, eedf, this->inelasticMatrix);
                integrate_source<LinIntegrator>(grid, *collision, eedf, this->inelasticMatrix);

                if (collision->isReverse())
                {
                    integrate_sup_sink<LinIntegrator>(grid, *collision, eedf, this->inelasticMatrix);
                    integrate_sup_source<LinIntegrator>(grid, *collision, eedf, this->inelasticMatrix);
                    // collision_integral_sup_sink(grid, eedf, this->inelasticMatrix, *collision);
                    // collision_integral_sup_source(grid, eedf, this->inelasticMatrix, *collision);
                }
            }
        }
    }

    // NOTE: This line is only required when the drift diffusion terms are also
    // divided by the cell width.
    // inelasticMatrix.array().rowwise() /= grid.duCells().array().transpose();
}

void equal_sharing_source(const Grid &grid, const Vector &eedf, Matrix &matrix, const EedfCollision &col)
{
    const double threshold = col.crossSection->threshold();
    const double target_density = col.getTarget()->delta();

    double u_start = grid.getNode(0) + threshold;

    // The grid iterator iterates over the source cells of the grid.
    Grid::Index i_grid = 0;

    // Find the target cell in which the left face of the first source cell lands.
    Grid::Index j_grid =
        std::upper_bound(grid.getNodes().begin(), grid.getNodes().end(), u_start) - grid.getNodes().begin() - 1;

    // Compute initial f_slope.
    double f_slope = u_start >= grid.getCell(j_grid) || j_grid == 0
                         ? std::log(eedf[j_grid + 1] / eedf[j_grid]) / grid.duNode(j_grid + 1)
                         : std::log(eedf[j_grid] / eedf[j_grid - 1]) / grid.duNode(j_grid);
    double u_f_start = u_start >= grid.getCell(j_grid) || j_grid == 0 ? grid.getCell(j_grid) : grid.getCell(j_grid - 1);

    // Get the raw cross section data.
    const auto &cs = col.crossSection->lookupTable().y();
    const auto &cs_energy = col.crossSection->lookupTable().x();

    Grid::Index j_cs = std::upper_bound(cs_energy.begin(), cs_energy.end(), u_start) - cs_energy.begin() - 1;

    double sig_slope = (cs[j_cs + 1] - cs[j_cs]) / (cs_energy[j_cs + 1] - cs_energy[j_cs]);
    double u_sig_start = cs_energy[j_cs];
    double u_sig_next = cs_energy[j_cs + 1];
    double sig_start = cs[j_cs];

    while (u_start < grid.uMax())
    {
        // If the current integration domain is outside the current cross
        // section cell, move to the next cell.
        if (u_start >= u_sig_next)
        {
            j_cs++;

            if (j_cs == cs.size() - 1)
            {
                sig_slope = 0.;
                u_sig_start = 0.;
                u_sig_next = std::numeric_limits<double>::max();
                sig_start = 0.;
            }
            else
            {
                sig_slope = (cs[j_cs + 1] - cs[j_cs]) / (cs_energy[j_cs + 1] - cs_energy[j_cs]);
                u_sig_start = cs_energy[j_cs];
                u_sig_next = cs_energy[j_cs + 1];
                sig_start = cs[j_cs];
            }
        }
        // If the current integration domain is outside the current source cell,
        // move to the next grid cell.
        if (u_start >= 2.0 * grid.getNode(i_grid + 1) + threshold)
        {
            i_grid++;
        }
        // If the current integration domain is outside the current target cell,
        // move to the next grid cell.
        if (u_start >= grid.getNode(j_grid + 1))
        {
            j_grid++;
            u_f_start = grid.getCell(j_grid);
        }
        // If the current integration domain starts at the current target cell
        // center, recompute the slope of the eedf.
        if (u_start == grid.getCell(j_grid) && j_grid < grid.nCells() - 1)
        {
            f_slope = std::log(eedf[j_grid + 1] / eedf[j_grid]) / grid.duNode(j_grid + 1);
        }
        // The end of the current integration domain is either the next source
        // cell center, the next cross section entry, the next target cell face,
        // or the grid boundary.
        const double u_end =
            std::min(std::min(u_sig_next, 2.0 * grid.getNode(i_grid + 1) + threshold),
                     u_start >= grid.getCell(j_grid) ? grid.getNode(j_grid + 1) : grid.getCell(j_grid));

        double integral_start = collision_integral(u_sig_start, sig_start, sig_slope, u_f_start, f_slope, u_start);
        double integral_end = collision_integral(u_sig_start, sig_start, sig_slope, u_f_start, f_slope, u_end);

        matrix(i_grid, j_grid) += 4.0 * target_density * (integral_end - integral_start);

        u_start = u_end;
    }
}

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

            // collision_integral_sink(grid, eedf, ionizationMatrix, *collision);
            integrate_sink<LinIntegrator>(grid, *collision, eedf, ionizationMatrix);

            switch (operatorType)
            {
            case IonizationOperatorType::conservative:
                // collision_integral_source(grid, eedf, ionizationMatrix, *collision);
                integrate_source<LinIntegrator>(grid, *collision, eedf, ionizationMatrix);
                break;
            case IonizationOperatorType::equalSharing:
                equal_sharing_source(grid, eedf, ionizationMatrix, *collision);
                break;
            }
        }
    }
}

SpatialGrowthOperator::SpatialGrowthOperator(const Grid &grid) : DriftDiffusionOperator(grid)
{
}

double mobility_integral(double u, double u_sig_start, double sig_start, double sig_slope, double u_f_start,
                         double u_f_end, double f_start, double f_div, double du)
{
    double value = 0.;
    // Linear
    if (sig_slope == 0.)
    {
        value = f_div * (1. / 6. * u * u / sig_start);
    }
    else
    {
        value = f_div / (3. * sig_slope * sig_slope) *
                (sig_slope * u +
                 (-sig_start + sig_slope * u_sig_start) * std::log(sig_start + sig_slope * (u - u_sig_start)));
    }
    // Logarithmic
    // NOTE: This is unstable.
    // if (sig_slope == 0.)
    // {
    //     value = std::exp((u - u_f_start) * f_div) * f_start * (f_div * u - 1.0) / (3 * f_div * sig_start);
    // }
    // else
    // {
    //     value = std::exp(-f_div * (sig_start / sig_slope + u_f_start)) * f_start / (3. * sig_slope * sig_slope) *
    //             (std::exp(f_div * (sig_start / sig_slope + u)) * sig_slope -
    //              std::exp(f_div * u_sig_start) * f_div * (sig_start - sig_slope * u_sig_start) *
    //                  std::expint(f_div * (sig_start / sig_slope + u - u_sig_start)));
    // }

    return value;
}

double integrate_mobility(const Grid &grid, const Vector &eedf, const Vector &total_cs)
{
    // The grid iterator iterates over the cells of the grid.
    Grid::Index i_grid = 0;

    // Compute initial f_slope using forward differences.
    double f_slope = (eedf[1] - eedf[0]) / (grid.getCell(1) - grid.getCell(0));
    // double f_slope = (std::log(eedf[1]) - std::log(eedf[0])) / (grid.getCell(1) - grid.getCell(0));
    double f_start = eedf[0];
    double u_f_start = grid.getCell(0);
    double u_f_end = grid.getCell(1);

    double sig_slope = (total_cs[1] - total_cs[0]) / (grid.getNode(1) - grid.getNode(0));
    double u_sig_start = 0.;
    double sig_start = total_cs[0];

    double u_start = 0;

    double mobility = 0;

    while (u_start < grid.uMax())
    {
        // If the current integration domain is outside the current grid cell,
        // move to the next grid cell.
        if (u_start >= grid.getNode(i_grid + 1))
        {
            i_grid++;
            u_f_start = grid.getCell(i_grid);
            u_f_end = grid.getCell(i_grid + 1);

            sig_slope = (total_cs[i_grid + 1] - total_cs[i_grid]) / (grid.getNode(i_grid + 1) - grid.getNode(i_grid));
            u_sig_start = grid.getNode(i_grid);
            sig_start = total_cs[i_grid];
        }
        // If the current integration domain starts at the current cell center,
        // recompute the slope of the eedf.
        if (u_start == grid.getCell(i_grid) && i_grid < grid.nCells() - 1)
        {
            f_slope = (eedf[i_grid + 1] - eedf[i_grid]) / (grid.getCell(i_grid + 1) - grid.getCell(i_grid));
            // f_slope = (std::log(eedf[i_grid + 1]) - std::log(eedf[i_grid])) /
            //           (grid.getCell(i_grid + 1) - grid.getCell(i_grid));
        }

        // The end of the current integration domain is either the next cell
        // center, the next cross section entry, or the grid boundary.
        const double u_end =
            std::min(u_start >= grid.getCell(i_grid) ? grid.getNode(i_grid + 1) : grid.getCell(i_grid), grid.uMax());

        mobility += mobility_integral(u_end, u_sig_start, sig_start, sig_slope, u_f_start, u_f_end, f_start, f_slope,
                                      u_end - u_start) -
                    mobility_integral(u_start, u_sig_start, sig_start, sig_slope, u_f_start, u_f_end, f_start, f_slope,
                                      u_end - u_start);

        u_start = u_end;
    }

    return -mobility;
}

void SpatialGrowthOperator::evaluate(const Grid &grid, const Vector &eedf, const Vector &total_cs, double EoN)
{
    Log<Message>::Warning("New: ", SI::gamma * integrate_mobility(grid, eedf, total_cs));
    const Vector Dnodes = grid.getNodes().array() / (3. * total_cs.array());
    const Vector D = (Dnodes.head(Dnodes.size() - 1) + Dnodes.tail(Dnodes.size() - 1)) / 2.0;
    Log<Message>::Warning("Old: ", -SI::gamma * fNodegPrimeEnergyIntegral(grid, Dnodes, eedf));
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
} // namespace experimental
} // namespace loki
