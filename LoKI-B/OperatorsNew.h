#include "LoKI-B/Grid.h"
#include "LoKI-B/LinearAlgebra.h"

namespace loki
{
namespace experimental
{
class DriftDiffusionOperator
{
  public:
    DriftDiffusionOperator(const Grid &grid);
    const Vector &drift_coefficient();
    const Vector &diffusion_coefficient();

  protected:
    Vector drift_coeff;
    Vector diff_coeff;
};
class ElasticOperator : public DriftDiffusionOperator
{
  public:
    ElasticOperator(const Grid &grid);
    void evaluate(const Grid &grid, const Vector &elasticCrossSection, double T_gas);
};

class FieldOperator : public DriftDiffusionOperator
{
  public:
    FieldOperator(const Grid &grid);
    void evaluate(const Grid &grid, const Vector &total_cs, double EoN);
};
} // namespace exp
} // namespace loki
