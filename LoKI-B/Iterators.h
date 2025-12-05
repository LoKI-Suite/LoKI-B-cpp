#ifndef ITERATORS_H
#define ITERATORS_H

#include "LoKI-B/LinearAlgebra.h"
#include <cmath>

namespace loki
{
// Forward declarations
class Grid;

class InterpolatingIterator
{
  public:
    InterpolatingIterator(const Vector &x, const Vector &y, double x_offset);

    bool shouldAdvance(double x_current) const;
    void advance();

    double linSlope() const;
    double logSlope() const;

    size_t index() const;
    double switchOn() const;
    double xLow() const;
    double xHigh() const;
    double yLow() const;
    double yHigh() const;

  private:
    const Vector &x_;
    const Vector &y_;
    Vector::Index current_index_;
    double switch_on_;
    double x_offset_;
    double x_low_, x_high_;
    double y_low_, y_high_;
};

class GridIterator
{
  public:
    GridIterator(const Grid &grid, double u_start);

    bool shouldAdvance(double x_current) const;
    void advance();

    size_t index() const;
    double xCell() const;
    double xLow() const;
    double xHigh() const;

  private:
    const Grid &grid_;
    size_t cell_index_;
};
} // namespace loki

#endif // ITERATORS_H
