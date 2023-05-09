#ifndef LOKI_CPP_GNUPLOT_H
#define LOKI_CPP_GNUPLOT_H

#include <string>
#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/Exports.h"

namespace loki {

lokib_export void writeGnuplot(std::ostream& os, const std::string &title, const std::string &xlabel, const std::string &ylabel, const Vector &x,
          const Vector &y);

} // namespace loki

#endif // LOKI_CPP_GNUPLOT_H
