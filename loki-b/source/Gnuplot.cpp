#include "LoKI-B/Gnuplot.h"

namespace loki {

void writeGnuplot(std::ostream& os, const std::string &title, const std::string &xlabel, const std::string &ylabel, const Vector &x,
          const Vector &y)
{
    os << "unset key" << std::endl;
    os << "set xlabel \"" << xlabel << "\"" << std::endl;
    os << "set ylabel \"" << ylabel << "\"" << std::endl;
    os << "set title \"" << title << "\"" << std::endl;
    os << "set xrange [" << x[0] << ":" << x[x.size() - 1] << "]" << std::endl;
    os << "set logscale y" << std::endl;
    // os << "set format y '%g'" << std::endl;
    os << "plot '-' w l" << std::endl;
    for (uint32_t i = 0; i < x.size(); ++i)
    {
        os << x[i] << "\t" << y[i] << '\n';
    }
    os << "e" << std::endl;
}

} // namespace loki
