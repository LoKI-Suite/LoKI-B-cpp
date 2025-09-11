#include "LoKI-B/Iterators.h"
#include "LoKI-B/LinearAlgebra.h"
namespace loki
{

class ConstIntegrator final
{
  public:
    static void set_rows(double x_low, double x_high, double x_offset, const GridIterator &grid_iter,
                         const InterpolatingIterator &cs_iter, const InterpolatingIterator &eedf_iter, Matrix &matrix);

  private:
    [[nodiscard]] static double value_at(double x, double offset, const InterpolatingIterator &cs_iter);
};

class LinearIntegrator
{
  public:
    static void set_rows(double x_low, double x_high, double x_offset, const GridIterator &grid_iter,
                         const InterpolatingIterator &cs_iter, const InterpolatingIterator &eedf_iter, Matrix &matrix);

  private:
    [[nodiscard]] static double first_term_at(double energy, double offset, const InterpolatingIterator &cs,
                                              const InterpolatingIterator &eedf);

    [[nodiscard]] static double second_term_at(double energy, double offset, const InterpolatingIterator &cs,
                                               const InterpolatingIterator &eedf);
};

class LogIntegrator
{
  public:
    static void set_rows(double x_low, double x_high, double x_offset, const GridIterator &grid_iter,
                         const InterpolatingIterator &cs_iter, const InterpolatingIterator &eedf_iter, Matrix &matrix);

  private:
    [[nodiscard]] static double value_at(double energy, double offset, const InterpolatingIterator &cs,
                                         const InterpolatingIterator &eedf);
};
} // namespace loki
