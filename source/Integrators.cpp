#include "LoKI-B/Integrators.h"
#include "LoKI-B/Iterators.h"
#include <cassert>
#include <cmath>

namespace loki
{
void LogIntegrator::setRows(double x_low, double x_high, double x_offset, double pre_factor,
                            const GridIterator &grid_iter, const InterpolatingIterator &cs_iter,
                            const InterpolatingIterator &eedf_iter, Eigen::MatrixXd &matrix)
{
    const auto row = static_cast<Eigen::Index>(grid_iter.index());
    const auto col = static_cast<Eigen::Index>(eedf_iter.index());

    matrix(row, col) +=
        pre_factor * (valueAt(x_high, x_offset, cs_iter, eedf_iter) - valueAt(x_low, x_offset, cs_iter, eedf_iter));
}

double LogIntegrator::valueAt(double energy, double offset, const InterpolatingIterator &cs,
                              const InterpolatingIterator &eedf)
{
    const double cs_slope = cs.linSlope();
    const double f_slope = eedf.logSlope();

    const double f_slope2 = f_slope * f_slope;
    const double f_slope3 = f_slope2 * f_slope;

    const double first = 2.0 * cs_slope / f_slope3;
    const double second = (energy + offset) * (cs.yLow() + (energy - cs.xLow()) * cs_slope) / f_slope;
    const double third = -(cs.yLow() + (2.0 * energy + offset - cs.xLow()) * cs_slope) / f_slope2;

    const double exponent = std::exp(f_slope * (energy - eedf.xLow()));

    return exponent * (first + second + third);
}

void LinIntegrator::setRows(double x_low, double x_high, double x_offset, double pre_factor,
                            const GridIterator &grid_iter, const InterpolatingIterator &cs_iter,
                            const InterpolatingIterator &eedf_iter, Eigen::MatrixXd &matrix)
{
    const auto row = static_cast<Eigen::Index>(grid_iter.index());
    const auto col0 = static_cast<Eigen::Index>(eedf_iter.index());
    const auto col1 = static_cast<Eigen::Index>(eedf_iter.index() + 1);

    matrix(row, col0) += pre_factor * (firstTermAt(x_high, x_offset, cs_iter, eedf_iter) -
                                       firstTermAt(x_low, x_offset, cs_iter, eedf_iter));
    matrix(row, col1) += pre_factor * (secondTermAt(x_high, x_offset, cs_iter, eedf_iter) -
                                       secondTermAt(x_low, x_offset, cs_iter, eedf_iter));
}

double LinIntegrator::firstTermAt(double energy, double offset, const InterpolatingIterator &cs,
                                  const InterpolatingIterator &eedf)
{
    const double cs_slope = cs.linSlope();

    const double e2 = energy * energy;
    const double e3 = e2 * energy;
    const double e4 = e2 * e2;

    const double first = -(e4 * cs_slope) / 4.0;
    const double second = -e3 * ((cs.yLow() - (eedf.xHigh() + cs.xLow()) * cs_slope) + offset * cs_slope) / 3.0;
    const double third = e2 *
                         (eedf.xHigh() * (cs.yLow() - cs.xLow() * cs_slope) -
                          offset * (cs.yLow() - (eedf.xHigh() + cs.xLow()) * cs_slope)) /
                         2.0;
    const double fourth = energy * eedf.xHigh() * offset * (cs.yLow() - cs.xLow() * cs_slope);

    const double denom = (eedf.xHigh() - eedf.xLow());
    assert(denom != 0.0 && "eedf.xHigh() - eedf.xLow() must be non-zero");
    return (first + second + third + fourth) / denom;
}

double LinIntegrator::secondTermAt(double energy, double offset, const InterpolatingIterator &cs,
                                   const InterpolatingIterator &eedf)
{
    const double cs_slope = cs.linSlope();

    const double e2 = energy * energy;
    const double e3 = e2 * energy;
    const double e4 = e2 * e2;

    const double first = e4 * cs_slope / 4.0;
    const double second = e3 * ((cs.yLow() - (eedf.xLow() + cs.xLow()) * cs_slope) + offset * cs_slope) / 3.0;
    const double third = -e2 *
                         (eedf.xLow() * (cs.yLow() - cs.xLow() * cs_slope) -
                          offset * (cs.yLow() - (eedf.xLow() + cs.xLow()) * cs_slope)) /
                         2.0;
    const double fourth = -energy * eedf.xLow() * offset * (cs.yLow() - cs.xLow() * cs_slope);

    const double denom = (eedf.xHigh() - eedf.xLow());
    assert(denom != 0.0 && "eedf.xHigh() - eedf.xLow() must be non-zero");
    return (first + second + third + fourth) / denom;
}

void ConstIntegrator::setRows(double x_low, double x_high, double x_offset, double pre_factor,
                              const GridIterator &grid_iter, const InterpolatingIterator &cs_iter,
                              const InterpolatingIterator &eedf_iter, Eigen::MatrixXd &matrix)
{
    const auto row = static_cast<Eigen::Index>(grid_iter.index());
    const auto col = static_cast<Eigen::Index>(eedf_iter.index());

    matrix(row, col) += pre_factor * (valueAt(x_high, x_offset, cs_iter) - valueAt(x_low, x_offset, cs_iter));
}

double ConstIntegrator::valueAt(double x, double offset, const InterpolatingIterator &cs_iter)
{
    const double cs_slope = cs_iter.linSlope();

    const double x2 = x * x;
    const double x3 = x2 * x;

    const double first = x3 * cs_slope / 3.0;
    const double second = x2 * (cs_iter.yLow() - (cs_iter.xLow() - offset) * cs_slope) / 2.0;
    const double third = x * offset * (cs_iter.yLow() - cs_iter.xLow() * cs_slope);

    return first + second + third;
}
} // namespace loki
