#include "LoKI-B/Enumeration.h"
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

class InelasticOperator
{
  public:
    InelasticOperator(const Grid &grid);
    void evaluate(const Grid &grid, const Vector &eedf, const EedfMixture &mixture);

    Matrix inelasticMatrix;
    Matrix superelasticMatrix;
};

class IonizationOperator
{
  public:
    IonizationOperator(const Grid &grid, IonizationOperatorType type);
    void evaluate(const Grid &grid, const Vector &eedf, const EedfMixture &mixture);

    const IonizationOperatorType operatorType;

    Matrix ionizationMatrix;
};

class SpatialGrowthOperator : public DriftDiffusionOperator
{
  public:
    SpatialGrowthOperator(const Grid &grid);
    void evaluate(const Grid &grid, const Vector &eedf, const Vector &total_cs, double EoN, const Matrix &ionizationMatrix, const Matrix &attachmentMatrix);
};

} // namespace experimental
} // namespace loki
