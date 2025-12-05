#include "LoKI-B/Iterators.h"
#include "LoKI-B/Grid.h"
#include <limits>

namespace loki
{
InterpolatingIterator::InterpolatingIterator(const Vector &x, const Vector &y, double x_offset)
    : x_(x), y_(y), current_index_(0), x_offset_(x_offset)
{
    switch_on_ = x_[1] - x_offset_;
    x_low_ = x_[0] - x_offset_;
    x_high_ = x_[1] - x_offset_;
    y_low_ = y_[0];
    y_high_ = y_[1];

    while (shouldAdvance(0.0))
    {
        advance();
    }
}

bool InterpolatingIterator::shouldAdvance(double x_current) const
{
    return x_current >= switch_on_;
}

void InterpolatingIterator::advance()
{
    if (current_index_ >= x_.size() - 2)
    {
        switch_on_ = std::numeric_limits<double>::max();
        return;
    }

    current_index_++;
    x_low_ = x_high_;
    y_low_ = y_high_;

    x_high_ = x_[current_index_ + 1] - x_offset_;
    y_high_ = y_[current_index_ + 1];

    switch_on_ = x_high_;
}

double InterpolatingIterator::linSlope() const
{
    return (y_high_ - y_low_) / (x_high_ - x_low_);
}

double InterpolatingIterator::logSlope() const
{
    return (std::log(y_high_) - std::log(y_low_)) / (x_high_ - x_low_);
}

size_t InterpolatingIterator::index() const
{
    return current_index_;
}
double InterpolatingIterator::switchOn() const
{
    return switch_on_;
}
double InterpolatingIterator::xLow() const
{
    return x_low_;
}
double InterpolatingIterator::xHigh() const
{
    return x_high_;
}
double InterpolatingIterator::yLow() const
{
    return y_low_;
}
double InterpolatingIterator::yHigh() const
{
    return y_high_;
}

GridIterator::GridIterator(const Grid &grid, double u_start) : grid_(grid), cell_index_(0)
{
    while (shouldAdvance(u_start))
    {
        advance();
    }
}

bool GridIterator::shouldAdvance(double x_current) const
{
    // NOTE: Removing the the addition of epsilon will cause the simulation to hang when integrating the superelastic
    // source term in some cases.
    return x_current + std::numeric_limits<double>::epsilon() >= grid_.getNode(cell_index_ + 1);
}

void GridIterator::advance()
{
    cell_index_++;
}

size_t GridIterator::index() const
{
    return cell_index_;
}
double GridIterator::xCell() const
{
    return grid_.getCell(cell_index_);
}
double GridIterator::xLow() const
{
    return grid_.getNode(cell_index_);
}
double GridIterator::xHigh() const
{
    return grid_.getNode(cell_index_ + 1);
}
} // namespace loki
