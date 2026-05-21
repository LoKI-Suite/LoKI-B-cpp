#include "LoKI-B/OperatorsNew.h"
#include "LoKI-B/EedfCollisions.h"
#include "LoKI-B/EedfMixture.h"
#include "LoKI-B/Enumeration.h"
#include "LoKI-B/Grid.h"
#include "LoKI-B/Integrators.h"
#include "LoKI-B/Iterators.h"
#include "LoKI-B/LinearAlgebra.h"
#include <cmath>
#include <stdexcept>

namespace loki
{
namespace experimental
{
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

// NOTE: this function implements the equal sharing source term using a variable substitution of `u' = 2*u + u_th`.
template <typename I>
void integrate_source_equal_sharing(const Grid &grid, const EedfCollision &col, const Vector &eedf, Matrix &mat)
{
    const auto target_density = col.getTarget()->delta();
    const auto &cs = *col.crossSection;

    const double u_offset = cs.threshold();
    // This is `u' = 2 * u + u_th`.
    double u_start = u_offset;

    GridIterator grid_iter(grid, 0.0);
    InterpolatingIterator cs_int(cs.lookupTable().x(), cs.lookupTable().y(), 0.0);
    InterpolatingIterator eedf_int(grid.getCells(), eedf, 0.0);

    while (u_start < grid.uMax())
    {
        while (grid_iter.shouldAdvance((u_start - u_offset) / 2.0))
        {
            grid_iter.advance();
        }
        while (cs_int.shouldAdvance(u_start))
        {
            cs_int.advance();
        }
        while (eedf_int.shouldAdvance(u_start))
        {
            eedf_int.advance();
        }

        double u_end =
            std::min({2. * grid_iter.xHigh() + u_offset, cs_int.switchOn(), eedf_int.switchOn(), grid.uMax()});

        // The additional division by `2` originates from the variable substitution.
        I::setRows(u_start, u_end, 0., 2. * target_density, grid_iter, cs_int, eedf_int, mat);

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
    superelasticMatrix.setZero();

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
                    integrate_sup_sink<LinIntegrator>(grid, *collision, eedf, this->superelasticMatrix);
                    integrate_sup_source<LinIntegrator>(grid, *collision, eedf, this->superelasticMatrix);
                }
            }
        }
    }

    // NOTE: This line is only required when the drift diffusion terms are also
    // divided by the cell width.
    inelasticMatrix.array().rowwise() /= grid.duCells().array().transpose();
}

IonizationOperator::IonizationOperator(const Grid &grid, IonizationOperatorType type)
    : operatorType(type), ionizationMatrix(grid.nCells(), grid.nCells()), includeNonConservativeIonization(false)
{
}

// NOTE: For now this only implements conservative ionization and equal sharing.
void IonizationOperator::evaluate(const Grid &grid, const Vector &eedf, const EedfMixture &mixture)
{
    bool has_valid_collisions = false;

    ionizationMatrix.setZero();

    for (const auto &cd : mixture.collision_data().data_per_gas())
    {
        for (const auto &collision : cd.collisions(CollisionType::ionization))
        {
            if (collision->crossSection->threshold() > grid.uMax())
                continue;

            has_valid_collisions = true;

            integrate_sink<LinIntegrator>(grid, *collision, eedf, ionizationMatrix);

            switch (operatorType)
            {
            case IonizationOperatorType::conservative:
                integrate_source<LinIntegrator>(grid, *collision, eedf, ionizationMatrix);
                break;
            case IonizationOperatorType::equalSharing:
                integrate_source_equal_sharing<LinIntegrator>(grid, *collision, eedf, ionizationMatrix);
                break;
            default:
                throw std::runtime_error(
                    "For now only conservative, and equal-sharing ionization modes are supported.");
            }
        }
    }
    ionizationMatrix.array().rowwise() /= grid.duCells().array().transpose();

    if (has_valid_collisions && operatorType != IonizationOperatorType::conservative)
        includeNonConservativeIonization = true;
}
} // namespace experimental
} // namespace loki
