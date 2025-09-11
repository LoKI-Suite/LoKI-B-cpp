#include "LoKI-B/Integrators.h"

namespace loki
{
double LogIntegrator::value_at(double energy, double offset, const InterpolatingIterator &cs,
                               const InterpolatingIterator &eedf)
{
    const double cs_slope = cs.lin_slope();
    const double f_slope = eedf.log_slope();

    const double inv_f = 1.0 / f_slope;
    const double inv_f2 = inv_f * inv_f;
    const double inv_f3 = inv_f2 * inv_f;

    const double first = 2.0 * cs_slope * inv_f3;
    const double second = (energy + offset) * (cs.y_low() + (energy - cs.x_low()) * cs_slope) * inv_f;
    const double third = -(cs.y_low() + (2.0 * energy + offset - cs.x_low()) * cs_slope) * inv_f2;

    return std::exp(f_slope * (energy - eedf.x_low())) * (first + second + third);
}

void LogIntegrator::set_rows(double x_low, double x_high, double x_offset, const GridIterator &grid_iter,
                             const InterpolatingIterator &cs_iter, const InterpolatingIterator &eedf_iter,
                             Matrix &matrix)
{
    matrix(grid_iter.index(), eedf_iter.index()) +=
        value_at(x_high, x_offset, cs_iter, eedf_iter) - value_at(x_low, x_offset, cs_iter, eedf_iter);
}

double ConstIntegrator::value_at(double energy, double offset, const InterpolatingIterator &cs)
{
    const double cs_slope = cs.lin_slope();

    const double x2 = energy * energy;
    const double x3 = x2 * energy;

    const double first = x3 * cs_slope / 3.0;
    const double second = x2 * (cs.y_low() - (cs.x_low() - offset) * cs_slope) / 2.0;
    const double third = energy * offset * (cs.y_low() - cs.x_low() * cs_slope);

    return first + second + third;
}

void ConstIntegrator::set_rows(double x_low, double x_high, double x_offset, const GridIterator &grid_iter,
                               const InterpolatingIterator &cs_iter, const InterpolatingIterator &eedf_iter,
                               Matrix &matrix)
{
    matrix(grid_iter.index(), eedf_iter.index()) +=
        value_at(x_high, x_offset, cs_iter) - value_at(x_low, x_offset, cs_iter);
}

double LinearIntegrator::first_term_at(double energy, double offset, const InterpolatingIterator &cs,
                                       const InterpolatingIterator &eedf)
{
    const double cs_slope = cs.lin_slope();

    const double e2 = energy * energy;
    const double e3 = e2 * energy;
    const double e4 = e2 * e2;

    const double first = -(e4 * cs_slope) / 4.0;
    const double second = -e3 * ((cs.y_low() - (eedf.x_high() + cs.x_low()) * cs_slope) + offset * cs_slope) / 3.0;
    const double third = e2 *
                         (eedf.x_high() * (cs.y_low() - cs.x_low() * cs_slope) -
                          offset * (cs.y_low() - (eedf.x_high() + cs.x_low()) * cs_slope)) /
                         2.0;
    const double fourth = energy * eedf.x_high() * offset * (cs.y_low() - cs.x_low() * cs_slope);

    return (first + second + third + fourth) / (eedf.x_high() - eedf.x_low());
}

double LinearIntegrator::second_term_at(double energy, double offset, const InterpolatingIterator &cs,
                                        const InterpolatingIterator &eedf)
{
    const double cs_slope = cs.lin_slope();

    const double e2 = energy * energy;
    const double e3 = e2 * energy;
    const double e4 = e2 * e2;

    const double first = e4 * cs_slope / 4.0;
    const double second = e3 * ((cs.y_low() - (eedf.x_low() + cs.x_low()) * cs_slope) + offset * cs_slope) / 3.0;
    const double third = -e2 *
                         (eedf.x_low() * (cs.y_low() - cs.x_low() * cs_slope) -
                          offset * (cs.y_low() - (eedf.x_low() + cs.x_low()) * cs_slope)) /
                         2.0;
    const double fourth = -energy * eedf.x_low() * offset * (cs.y_low() - cs.x_low() * cs_slope);

    return (first + second + third + fourth) / (eedf.x_high() - eedf.x_low());
}

void LinearIntegrator::set_rows(double x_low, double x_high, double x_offset, const GridIterator &grid_iter,
                                const InterpolatingIterator &cs_iter, const InterpolatingIterator &eedf_iter,
                                Matrix &matrix)
{
    const auto r = grid_iter.index();
    const auto c = eedf_iter.index();

    matrix(r, c) +=
        first_term_at(x_high, x_offset, cs_iter, eedf_iter) - first_term_at(x_low, x_offset, cs_iter, eedf_iter);
    matrix(r, c + 1) +=
        second_term_at(x_high, x_offset, cs_iter, eedf_iter) - second_term_at(x_low, x_offset, cs_iter, eedf_iter);
}
} // namespace loki
