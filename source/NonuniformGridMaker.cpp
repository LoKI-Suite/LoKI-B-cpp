#include "LoKI-B/NonuniformGridMaker.h"

namespace loki
{
loki::Vector writeGrid(const loki::Vector nodeDistribution, const double diff)
{
    double dumin = (nodeDistribution.tail(nodeDistribution.size() - 1) - nodeDistribution.head(nodeDistribution.size() - 1)).minCoeff();
    loki::Vector u(3000);
    double lowdiff = diff - 1.;
    u[0] = 0.0;
    u[1] = dumin;
    int n = 1;
    if ((nodeDistribution[1] - dumin)/dumin > diff)
    {
        while ((nodeDistribution[1] - u[n])/(u[n] - u[n - 1]) > diff)
        {
            u[n+1] = (nodeDistribution[1] + u[n])/2.;
            n++;
        }
    }
    for (loki::Grid::Index k = 2; k < nodeDistribution.size(); ++k)
    {
        if ((nodeDistribution[k]- u[n])/(u[n] - u[n-1]) > diff)
        {
            while ((nodeDistribution[k]-u[n])/(u[n] - u[n - 1]) > diff)
            {
                u[n + 1] = (nodeDistribution[k] + u[n])/2.;
                n++;

            }
        }

        if ((nodeDistribution[k]- u[n])/(u[n] - u[n-1]) < lowdiff)
        {
            continue;
        }
        u[n + 1] =  nodeDistribution[k];
        n++;
    }
    n++;

    u.conservativeResize(n);

    return u;
}

void removeDuplicates(std::vector<double> &flatEnergies)
{
    std::sort(flatEnergies.begin(), flatEnergies.end());
    auto it = std::unique(flatEnergies.begin(), flatEnergies.end());
    flatEnergies.erase(it, flatEnergies.end());
}

loki::Vector vectorToVector(std::vector<double> &flatEnergies)
{
    Eigen::Map<loki::Vector> flatGrid(flatEnergies.data(), flatEnergies.size());

    return flatGrid;
}

loki::Vector removeCloseNeighbors(loki::Vector &gridNodes, double alpha)
{
    int newSize = 1;
    loki::Vector filteredNodes(gridNodes.size());
    filteredNodes[0] = gridNodes[0];
    for (int i = 1; i < gridNodes.size(); ++i) {
        if (std::abs(gridNodes[i] - gridNodes[i - 1]) >= alpha) {
            filteredNodes[newSize++] = gridNodes[i]; 
        }
    }
    filteredNodes.conservativeResize(newSize);

    return filteredNodes;
}

loki::Vector makeGridFromMixture(const loki::EedfMixture &mixture, const double maxEnergy, const double closestNeighbors, const double cellgrowth)
{
    using namespace loki;
    std::vector<double> flatData;

    for (uint32_t i = 0; i < mixture.collision_data().m_collisions.size(); i++)
    {
        for (uint32_t j = 0; j < mixture.collision_data().m_collisions[i].m_coll->crossSection->lookupTable().x().size(); j++)
        {
            flatData.push_back(mixture.collision_data().m_collisions[i].m_coll->crossSection->lookupTable().x()[j]);
        }
    }
    removeDuplicates(flatData);
    Vector gridNodes = vectorToVector(flatData);
    Vector filteredGrid = removeCloseNeighbors(gridNodes, closestNeighbors);
    Vector grid = writeGrid(filteredGrid, cellgrowth);

    return grid;
}
}