#pragma once

#include <Eigen/Dense>

namespace loki
{
class InterpolatingIterator
{
  public:
    using index_t = Eigen::Index;

    /// Construct with references to Eigen::VectorXd (non-owning).
    /// Precondition: x.size() == y.size() and size >= 2.
    explicit InterpolatingIterator(const Eigen::VectorXd &x, const Eigen::VectorXd &y, double x_offset = 0.0);

    /// Advance as long as x_current >= switch_on().
    void advance_to(double x_current);

    /// Linear slope over the current segment.
    [[nodiscard]] double lin_slope() const noexcept;

    /// Log-linear slope over the current segment.
    [[nodiscard]] double log_slope() const;

    // Accessors 
    [[nodiscard]] index_t index() const noexcept;
    [[nodiscard]] double switch_on() const noexcept;

    [[nodiscard]] double x_low() const noexcept;
    [[nodiscard]] double x_high() const noexcept;
    [[nodiscard]] double y_low() const noexcept;
    [[nodiscard]] double y_high() const noexcept;

    [[nodiscard]] double x_offset() const noexcept;

  private:
    const Eigen::VectorXd &m_x;
    const Eigen::VectorXd &m_y;

    index_t m_current_index;
    double m_x_offset;
    double m_switch_on;

    double m_x_low;
    double m_x_high;
    double m_y_low;
    double m_y_high;
};
} // namespace loki
