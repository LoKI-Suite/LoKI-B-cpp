#include "LoKI-B/Iterators.h"
#include "LoKI-B/Log.h"

#include <cassert>
#include <cmath>
#include <limits>

namespace loki
{
InterpolatingIterator::InterpolatingIterator(const Eigen::VectorXd &x, const Eigen::VectorXd &y, double x_offset)
    : m_x(x), m_y(y), m_current_index(0), m_x_offset(x_offset), m_switch_on(0.0), m_x_low(0.0), m_x_high(0.0),
      m_y_low(0.0), m_y_high(0.0)
{
    if (m_x.size() != m_y.size())
    {
        Log<Message>::Error("InterpolatingIterator: x and y must have the same length.");
    }
    if (m_x.size() < 2)
    {
        Log<Message>::Error("InterpolatingIterator: x and y must contain at least two elements.");
    }

    m_x_low = m_x(0) - m_x_offset;
    m_x_high = m_x(1) - m_x_offset;
    m_y_low = m_y(0);
    m_y_high = m_y(1);
    m_switch_on = m_x_high;

    advance_to(0.0);
}

void InterpolatingIterator::advance_to(double x_current)
{
    const auto n = m_x.size();

    while (x_current >= m_switch_on)
    {
        if (m_current_index >= n - 2)
        {
            // No more segments; keep switch_on at infinity.
            m_switch_on = std::numeric_limits<double>::infinity();
            break;
        }

        ++m_current_index;

        m_x_low = m_x_high;
        m_y_low = m_y_high;
        m_x_high = m_x(m_current_index + 1) - m_x_offset;
        m_y_high = m_y(m_current_index + 1);

        m_switch_on = m_x_high;
    }
}

double InterpolatingIterator::lin_slope() const noexcept
{
    assert((m_x_high - m_x_low) != 0.0 && "Zero-length x segment: slope undefined.");
    return (m_y_high - m_y_low) / (m_x_high - m_x_low);
}

double InterpolatingIterator::log_slope() const
{
    assert(m_y_low > 0.0 && "log_slope requires y_low > 0");
    assert(m_y_high > 0.0 && "log_slope requires y_high > 0");
    assert((m_x_high - m_x_low) != 0.0 && "Zero-length x segment: slope undefined.");
    return (std::log(m_y_high) - std::log(m_y_low)) / (m_x_high - m_x_low);
}

InterpolatingIterator::index_t InterpolatingIterator::index() const noexcept
{
    return m_current_index;
}

double InterpolatingIterator::switch_on() const noexcept
{
    return m_switch_on;
}

double InterpolatingIterator::x_low() const noexcept
{
    return m_x_low;
}

double InterpolatingIterator::x_high() const noexcept
{
    return m_x_high;
}

double InterpolatingIterator::y_low() const noexcept
{
    return m_y_low;
}

double InterpolatingIterator::y_high() const noexcept
{
    return m_y_high;
}

double InterpolatingIterator::x_offset() const noexcept
{
    return m_x_offset;
}
} // namespace loki
