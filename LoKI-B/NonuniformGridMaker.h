#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/EedfMixture.h"

namespace loki
{
    loki::Vector writeGrid(const loki::Vector nodeDistribution, const double diff);

    void removeDuplicates(std::vector<double> &flatEnergies);

    loki::Vector vectorToVector(std::vector<double> &flatEnergies);

    loki::Vector removeCloseNeighbors(loki::Vector &gridNodes, double alpha);

    loki::Vector makeGridFromMixture(const loki::EedfMixture &mixture, const double maxEnergy, const double closestNeighbors, const double cellgrowth);
}