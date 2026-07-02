#ifndef LOKIB_CPP_OPERATORSNEW_H
#define LOKIB_CPP_OPERATORSNEW_H

#include "LoKI-B/Enumeration.h"
#include "LoKI-B/Grid.h"
#include "LoKI-B/LinearAlgebra.h"

namespace loki
{
namespace experimental
{
template<typename Integrator>
class InelasticOperator
{
  public:
    InelasticOperator(const Grid &grid);
    void evaluate(const Grid &grid, const Vector &eedf, const EedfMixture &mixture);

    Matrix inelasticMatrix;
    Matrix superelasticMatrix;
};

template<typename Integrator>
class IonizationOperator
{
  public:
    IonizationOperator(const Grid &grid, IonizationOperatorType type);
    void evaluate(const Grid &grid, const Vector &eedf, const EedfMixture &mixture);

    const IonizationOperatorType operatorType;

    Matrix ionizationMatrix;
    /* NOTE: the following is not a configuration parameter, but the
     * result of introspection of the reaction lists. The results also depend
     * on uMax.
     */
    bool includeNonConservativeIonization;
};
} // namespace experimental
} // namespace loki

#endif // LOKIB_CPP_ITERATORS_H
