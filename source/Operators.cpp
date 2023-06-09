#include "LoKI-B/Operators.h"
#include "LoKI-B/Constant.h"

#include <cassert>

namespace loki {

CAROperator::CAROperator(const CARGases& cg)
: carGases(cg)
{
}

void CAROperator::evaluate(const Grid& grid)
{
    const double a02 = Constant::bohrRadius*Constant::bohrRadius;
    double sigma0B = 0.;
    for (const auto &gas : carGases)
    {
        const double Qau = gas->electricQuadrupoleMoment/(Constant::electronCharge*a02);
        sigma0B += gas->fraction * Qau * Qau * gas->rotationalConstant;
    }
    sigma0B *= (8.*Constant::pi*a02/15.)*sigma0B;

    g = grid.getNodes() * (4. * sigma0B);
    // Boundary conditions. See Tejero2019 below equation 16b.
    g[0] = 0.;
    g[grid.nCells()] = 0.;
}

void CAROperator::evaluate(const Grid& grid, double Tg, SparseMatrix& mat)
{
    // update g
    evaluate(grid);

    const double c_CAR = Constant::kBeV * Tg;
    const double factor1 = (c_CAR / grid.du() + 0.5) / grid.du();
    const double factor2 = (c_CAR / grid.du() - 0.5) / grid.du();
    for (Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        mat.coeffRef(k, k) = -(g[k] * factor1 + g[k + 1] * factor2);

        if (k > 0)
            mat.coeffRef(k, k - 1) = g[k] * factor2;

        if (k < grid.nCells() - 1)
            mat.coeffRef(k, k + 1) = g[k + 1] * factor1;
    }
}

void CAROperator::evaluatePower(const Grid& grid, const Vector& eedf, double Tg, double& net, double& gain, double& loss) const
{
    /** \todo See what is the best way to write this. Note that
     *    (g[k + 1] * auxLow - g[k] * auxHigh)
     *  = (g[k + 1] - g[k]) * kTg - (g[k + 1] + g[k]) * du/2,
     *  so
     *  net = gain - sum_k eedf[k] (g[k + 1] + g[k]) * du/2 := gain + loss
     */
    const double kTg = Constant::kBeV * Tg;
    const double auxHigh = kTg + grid.du() * .5; // aux1
    const double auxLow  = kTg - grid.du() * .5; // aux2
    net = 0.0;
    gain = 0.0;
    for (Grid::Index k = 0; k < grid.nCells() - 1; ++k)
    {
        net += eedf[k] * (g[k + 1] * auxLow - g[k] * auxHigh);
        gain += eedf[k] * (g[k + 1] - g[k]);
    }
    net *= SI::gamma;
    gain *= SI::gamma * kTg;
    loss = net - gain;
}

ElasticOperator::ElasticOperator()
{
}

void ElasticOperator::evaluate(const Grid& grid, const Vector& elasticCrossSection)
{
    g = grid.getNodes().cwiseAbs2().cwiseProduct(elasticCrossSection) * 2;
    g[0] = 0.;
    g[g.size() - 1] = 0.;
}

void ElasticOperator::evaluate(const Grid& grid, const Vector& elasticCrossSection, double Tg, SparseMatrix& mat)
{
    // update g
    evaluate(grid,elasticCrossSection);

    const double c_el = Constant::kBeV * Tg;

    const double factor1 = (c_el / grid.du() + 0.5) / grid.du();
    const double factor2 = (c_el / grid.du() - 0.5) / grid.du();
    for (Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        mat.coeffRef(k, k) = -(g[k] * factor1 + g[k + 1] * factor2);

        if (k > 0)
            mat.coeffRef(k, k - 1) = g[k] * factor2;

        if (k < grid.nCells() - 1)
            mat.coeffRef(k, k + 1) = g[k + 1] * factor1;
    }
}

void ElasticOperator::evaluatePower(const Grid& grid, const Vector& eedf, double Tg, double& net, double& gain, double& loss) const
{
    const double kTg = Constant::kBeV * Tg;
    const double auxHigh = kTg + grid.du() * .5; // aux1
    const double auxLow  = kTg - grid.du() * .5; // aux2
    net = 0.0;
    gain = 0.0;
    for (Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        net += eedf[k] * (g[k + 1] * auxLow - g[k] * auxHigh);
        gain += eedf[k] * (g[k + 1] - g[k]);
    }
    net *= SI::gamma;
    gain *= SI::gamma * kTg;
    loss = net-gain;
}

FieldOperator::FieldOperator(const Grid& grid)
 : g(grid.getNodes().size())
{
}

void FieldOperator::evaluate(const Grid& grid, const Vector& totalCS, double EoN, double WoN)
{
    assert(g.size()==grid.getNodes().size());
    const double me = Constant::electronMass;
    const double e = Constant::electronCharge;
    g[0] = 0.;
    for (Grid::Index i=1; i!= g.size()-1; ++i)
    {
        g[i] = (EoN * EoN / 3) * grid.getNode(i) /
          (totalCS[i] + (me * WoN * WoN / (2 * e)) / (grid.getNode(i)*totalCS[i]));
    }
    g[g.size() - 1] = 0.;
}

void FieldOperator::evaluate(const Grid& grid, const Vector& totalCS, double EoN, double WoN, SparseMatrix& mat)
{
    // update g
    evaluate(grid,totalCS,EoN,WoN);
    const double sqStep = grid.du() * grid.du();

    for (Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        mat.coeffRef(k, k) = -(g[k] + g[k + 1]) / sqStep;

        if (k > 0)
            mat.coeffRef(k, k - 1) = g[k] / sqStep;

        if (k < grid.nCells() - 1)
            mat.coeffRef(k, k + 1) = g[k + 1] / sqStep;
    }
}

void FieldOperator::evaluatePower(const Grid& grid, const Vector& eedf, double& power) const
{
    power = 0;
    for (Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        power += eedf[k] * (g[k + 1] - g[k]);
    }
    power *= SI::gamma;
}

ElectronElectronOperator::ElectronElectronOperator()
{
    alphaEE=0.0;
}

void ElectronElectronOperator::initialize(const Grid& grid)
{
    A.setZero(grid.nCells());
    B.setZero(grid.nCells());
}

void ElectronElectronOperator::evaluatePower(const Grid& grid, const Vector& eedf, double& power) const
{
        power = (-SI::gamma * grid.du() * grid.du()) * (A - B).dot(eedf);
        //std::cout << "EE POWER: " << power << std::endl;
}

} // namespace loki
