
#ifndef INTEGRATORS_H
#define INTEGRATORS_H

#include "LoKI-B/LinearAlgebra.h"

namespace loki
{
// Forward declarations
class GridIterator;
class InterpolatingIterator;

class LogIntegrator
{
  public:
    static void setRows(double x_low, double x_high, double x_offset, double pre_factor, const GridIterator &grid_iter,
                        const InterpolatingIterator &cs_iter, const InterpolatingIterator &eedf_iter,
                        loki::Matrix &matrix);

  private:
    static double valueAt(double energy, double offset, const InterpolatingIterator &cs,
                          const InterpolatingIterator &eedf);
};

class LinIntegrator
{
  public:
    static void setRows(double x_low, double x_high, double x_offset, double pre_factor, const GridIterator &grid_iter,
                        const InterpolatingIterator &cs_iter, const InterpolatingIterator &eedf_iter,
                        loki::Matrix &matrix);

  private:
    static double firstTermAt(double energy, double offset, const InterpolatingIterator &cs,
                              const InterpolatingIterator &eedf);

    static double secondTermAt(double energy, double offset, const InterpolatingIterator &cs,
                               const InterpolatingIterator &eedf);
};

class ConstIntegrator
{
  public:
    static void setRows(double x_low, double x_high, double x_offset, double pre_factor, const GridIterator &grid_iter,
                        const InterpolatingIterator &cs_iter, const InterpolatingIterator &eedf_iter,
                        loki::Matrix &matrix);

  private:
    static double valueAt(double x, double offset, const InterpolatingIterator &cs_iter);
};
} // namespace loki

#endif // INTEGRATORS_H
