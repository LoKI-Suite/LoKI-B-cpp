#ifndef LOKI_CPP_IDEAS_HESSENBERG_EXTRA_H
#define LOKI_CPP_IDEAS_HESSENBERG_EXTRA_H

#include "LoKI-B/LinearAlgebra.h"

namespace loki {
namespace LinAlg {

/** \todo Document these functions, add tests.
 */

void hessenbergReductionPartialPiv(double *A, const uint32_t *c, uint32_t *p, uint32_t n, uint32_t cn);

double *hessenberg(double *A, double *b, const uint32_t *p, uint32_t n, HessenbergWorkspace& hws);

double *hessenberg(double *A, double *b, const uint32_t *p, uint32_t n);

} // namespace LinAlg
} // namespace loki

#endif // LOKI_CPP_IDEAS_HESSENBERG_EXTRA_H
