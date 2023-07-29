#include "LoKI-B/EedfUtilities.h"
#include <cmath>

namespace loki {

void normalizeEDF(Vector& edf, const Grid& grid)
{
    edf /= edf.dot(grid.getCells().cwiseSqrt() * grid.du());
}

void makePrescribedEDF(Vector& edf, const Grid& grid, double g, double T_eV)
{
    edf.resize(grid.nCells());
    const double gamma_3_2g = std::tgamma(3/(2*g));
    const double gamma_5_2g = std::tgamma(5/(2*g));
    /** \todo This follows the MATLAB code. But the prefactors can be removed,
     *  since this is normalized afterwards anyway.
     *  \todo The risk of underflows should be investigated. (And diagnosed?)
     */
    for (Grid::Index k = 0; k != grid.nCells(); ++k)
    {
        const double energy = grid.getCell(k);
        edf(k) = std::pow(gamma_5_2g,1.5) * std::pow(gamma_3_2g,-2.5) * std::pow(2/(3*T_eV),1.5)
            *std::exp(-std::pow(energy*gamma_5_2g*2/(gamma_3_2g*3*T_eV),g));
    }
    normalizeEDF(edf,grid);
}

Vector makePrescribedEDF(const Grid& grid, double g, double T_eV)
{
    Vector eedf(grid.nCells());
    makePrescribedEDF(eedf,grid,g,T_eV);
    return eedf;
}

} // namespace loki

