#include "LoKI-B/OperatorsNew.h"
#include "LoKI-B/Constant.h"
#include "LoKI-B/Grid.h"

namespace loki
{
namespace experimental
{
DriftDiffusionOperator::DriftDiffusionOperator(const Grid &grid)
    : drift_coeff(grid.getNodes().size()), diff_coeff(grid.getNodes().size())
{
}
const Vector &DriftDiffusionOperator::drift_coefficient()
{
    return this->drift_coeff;
}
const Vector &DriftDiffusionOperator::diffusion_coefficient()
{
    return this->diff_coeff;
}

ElasticOperator::ElasticOperator(const Grid &grid) : DriftDiffusionOperator(grid)
{
}
void ElasticOperator::evaluate(const Grid &grid, const Vector &elasticCrossSection, double T_gas)
{
    this->drift_coeff = grid.getNodes().cwiseAbs2().cwiseProduct(elasticCrossSection);
    this->diff_coeff = -(Constant::kBeV * T_gas) * this->drift_coefficient();
}

FieldOperator::FieldOperator(const Grid &grid) : DriftDiffusionOperator(grid)
{
}
void FieldOperator::evaluate(const Grid &grid, const Vector &total_cs, double EoN)
{
    this->drift_coeff.setZero();
    this->diff_coeff = -1. / 3. * (EoN * EoN) * grid.getNodes().cwiseQuotient(total_cs);
}
} // namespace exp
} // namespace loki
