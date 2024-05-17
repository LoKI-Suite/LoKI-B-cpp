#include "LoKI-B/Operators.h"
#include "LoKI-B/Constant.h"

#include <cassert>

namespace loki {

CAROperator::CAROperator(const CARGases& cg)
: carGases(cg), m_sigma0B(0)
{
    /** \todo Check here also that electricQuadrupoleMoment
     *  and rotationalConstant have been set for all CAR gases?
     */
}

void CAROperator::evaluate(const Grid& grid)
{
    constexpr const double a02 = Constant::bohrRadius*Constant::bohrRadius;
    m_sigma0B = 0.;
    for (const auto &gas : carGases)
    {
        const double Qau = gas->electricQuadrupoleMoment/(Constant::electronCharge*a02);
        m_sigma0B += gas->fraction * Qau * Qau * gas->rotationalConstant;
    }
    m_sigma0B *= (8.*Constant::pi*a02/15.);

    g = grid.getNodes() * (4. * m_sigma0B);
    // Boundary conditions. See Tejero2019 below equation 16b.
    g[0] = 0.;
    g[grid.nCells()] = 0.;
}

void CAROperator::evaluate(const Grid& grid, double Tg, SparseMatrix& mat)
{
    // update g
    evaluate(grid);
    const double c_CAR = Constant::kBeV * Tg;

    if (grid.isUniform())
    {
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

    } else
    {   
        for (Grid::Index k = 0; k < grid.nCells(); ++k)
        {
            mat.coeffRef(k, k) = 0;

            if (k > 0)
            {   
                const double Bmin = -grid.duCell(k) / grid.duNode(k) / 2 + c_CAR/grid.duNode(k);
                const double Amin = grid.duCell(k-1) / grid.duNode(k) / 2 + c_CAR/grid.duNode(k);
                
                mat.coeffRef(k, k - 1) = g[k] * Bmin / grid.duCell(k);
                mat.coeffRef(k, k) += -g[k] * Amin / grid.duCell(k);
            }
            
            if (k < grid.nCells() - 1)
            {
                const double Bplus = -grid.duCell(k+1) / grid.duNode(k+1) / 2 + c_CAR/grid.duNode(k+1);
                const double Aplus = grid.duCell(k) / grid.duNode(k+1) / 2 + c_CAR/grid.duNode(k+1);

                mat.coeffRef(k, k + 1) = g[k + 1] * Aplus / grid.duCell(k);
                mat.coeffRef(k, k) += -g[k + 1] * Bplus / grid.duCell(k);
            }
        }
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
    if (grid.isUniform())
    {
        double auxHigh = kTg + grid.du() * .5; // aux1
        double auxLow  = kTg - grid.du() * .5; // aux2
        net = 0.0;
        gain = 0.0;
        for (Grid::Index k = 0; k < grid.nCells(); ++k)
        {
            net += eedf[k] * (g[k + 1] * auxLow - g[k] * auxHigh);
            gain += eedf[k] * (g[k + 1] - g[k]);
        }
    } else
    {
        const auto n = grid.nCells();

        Vector auxHigh = kTg * Vector::Ones(n + 1) + grid.duNodes() * .5; // aux1
        Vector auxLow  = kTg * Vector::Ones(n + 1) - grid.duNodes() * .5; // aux2

        net = (eedf.array() * (g.tail(n).array() * auxLow.tail(n).array() - g.head(n).array() * auxHigh.head(n).array())).sum();
        gain = (eedf.array() * (g.tail(n) - g.head(n)).array()).sum();
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
    
    if (grid.isUniform())
    {
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
    } else
    {   
        for (Grid::Index k = 0; k < grid.nCells(); ++k)
        {
            mat.coeffRef(k, k) = 0;

            if (k > 0)
            {   
                const double Bmin = -grid.duCell(k) / grid.duNode(k) / 2 + c_el/grid.duNode(k);
                const double Amin = grid.duCell(k-1) / grid.duNode(k) / 2 + c_el/grid.duNode(k);

                mat.coeffRef(k, k - 1) = g[k] * Bmin / grid.duCell(k);
                mat.coeffRef(k, k) += -g[k] * Amin / grid.duCell(k);
            }
            
            if (k < grid.nCells() - 1)
            {
                const double Bplus = -grid.duCell(k+1) / grid.duNode(k+1) / 2 + c_el/grid.duNode(k+1);
                const double Aplus = grid.duCell(k) / grid.duNode(k+1) / 2 + c_el/grid.duNode(k+1);

                mat.coeffRef(k, k + 1) = g[k + 1] * Aplus / grid.duCell(k);
                mat.coeffRef(k, k) += -g[k + 1] * Bplus / grid.duCell(k);
            }
        }
    }
}

void ElasticOperator::evaluatePower(const Grid& grid, const Vector& eedf, double Tg, double& net, double& gain, double& loss) const
{
    const double kTg = Constant::kBeV * Tg;
    if (grid.isUniform())
    {
        double auxHigh = kTg + grid.du() * .5; // aux1
        double auxLow  = kTg - grid.du() * .5; // aux2
        net = 0.0;
        gain = 0.0;
        for (Grid::Index k = 0; k < grid.nCells(); ++k)
        {
            net += eedf[k] * (g[k + 1] * auxLow - g[k] * auxHigh);
            gain += eedf[k] * (g[k + 1] - g[k]);
        }
    } else
    {
        const auto n = grid.nCells();

        Vector auxHigh = kTg * Vector::Ones(n + 1) + grid.duNodes() * .5; // aux1
        Vector auxLow  = kTg * Vector::Ones(n + 1) - grid.duNodes() * .5; // aux2

        net = (eedf.array() * (g.tail(n).array() * auxLow.tail(n).array() - g.head(n).array() * auxHigh.head(n).array())).sum();
        gain = (eedf.array() * (g.tail(n) - g.head(n)).array()).sum();
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
    g[0] = 0.;
    for (Grid::Index i=1; i!= g.size()-1; ++i)
    {
        g[i] = (EoN * EoN / 3) * grid.getNode(i) /
          (totalCS[i] + ( WoN * WoN / (SI::gamma*SI::gamma)) / (grid.getNode(i)*totalCS[i]));
    }
    g[g.size() - 1] = 0.;
}

void FieldOperator::evaluate(const Grid& grid, const Vector& totalCS, double EoN, double WoN, SparseMatrix& mat)
{
    // update g
    evaluate(grid,totalCS,EoN,WoN);

    if (grid.isUniform())
    {
        const double sqStep = grid.du() * grid.du();

        for (Grid::Index k = 0; k < grid.nCells(); ++k)
        {
            mat.coeffRef(k, k) = -(g[k] + g[k + 1]) / sqStep;

            if (k > 0)
                mat.coeffRef(k, k - 1) = g[k] / sqStep;

            if (k < grid.nCells() - 1)
                mat.coeffRef(k, k + 1) = g[k + 1] / sqStep;
        }
    } else
    {
        for (Grid::Index k = 0; k < grid.nCells(); ++k)
        {
            mat.coeffRef(k, k) = 0;

            if (k > 0)
            {
                const double Amin = 1/grid.duNode(k);
                const double Bmin = 1/grid.duNode(k);

                mat.coeffRef(k, k - 1) = g[k] * Bmin / grid.duCell(k);
                mat.coeffRef(k, k) += -g[k] * Amin / grid.duCell(k);
            }
            if (k < grid.nCells() - 1)
            {
                const double Aplus = 1/grid.duNode(k+1);
                const double Bplus = 1/grid.duNode(k+1);

                mat.coeffRef(k, k + 1) = g[k + 1] * Aplus / grid.duCell(k);
                mat.coeffRef(k, k) += -g[k + 1] * Bplus / grid.duCell(k);
            }
        }
    }
}

void FieldOperator::evaluatePower(const Grid& grid, const Vector& eedf, double& power) const
{
    const auto n = grid.nCells();
    power = SI::gamma * eedf.cwiseProduct(g.tail(n) - g.head(n)).sum();
}

std::array<double, 2> alphaDistribution(double targetCell, double uMin, double uPlus, double frac = 1.)
{
    std::array<double, 2> alpha;
    alpha[0] = frac*(uPlus - targetCell) / (uPlus - uMin);
    alpha[1] = frac - alpha[0];

    return alpha;
}

std::vector<std::tuple<int, double>> distributeOneCell(const Grid& grid, double targetCell, int targetBegin)
{
    std::vector<std::tuple<int, double>> alpha;
    std::array<double, 2> alphaMinPlus;

    if (targetCell > grid.getCell(targetBegin))
    {
        alphaMinPlus = alphaDistribution(targetCell, grid.getCell(targetBegin), grid.getCell(targetBegin + 1));
        alpha.push_back(std::make_tuple(targetBegin, alphaMinPlus[0]));
        alpha.push_back(std::make_tuple(targetBegin + 1, alphaMinPlus[1]));
    } else if (targetCell < grid.getCell(targetBegin))
    {
        alphaMinPlus = alphaDistribution(targetCell, grid.getCell(targetBegin - 1), grid.getCell(targetBegin));
        alpha.push_back(std::make_tuple(targetBegin - 1, alphaMinPlus[0]));
        alpha.push_back(std::make_tuple(targetBegin, alphaMinPlus[1]));
    } else
    {
        alpha.push_back(std::make_tuple(targetBegin, 1.));
    }

    return alpha;
}

std::vector<std::tuple<int, double>> distributeTwoCells(const Grid& grid, double targetCell, int targetBegin, int targetEnd)
{
    std::vector<std::tuple<int, double>> alpha;
    std::array<double, 2> alphaMinPlus = alphaDistribution(targetCell, grid.getCell(targetBegin), grid.getCell(targetEnd));
    alpha.push_back(std::make_tuple(targetBegin, alphaMinPlus[0]));
    alpha.push_back(std::make_tuple(targetEnd, alphaMinPlus[1]));
    return alpha;
}

std::vector<std::tuple<int, double>> distributeNCells(const Grid& grid, double targetCell, int targetBegin, int targetEnd, 
                                                     Grid::Index origin, double threshold, bool reverse = false, double frac = 1.0)
{
    double targetMiddleLeft;
    double targetMiddleRight;
    double fractionLeft;
    double fractionRight;
    if (reverse)
    {
        targetMiddleLeft = (grid.getNode(targetBegin + 1) + frac * grid.getNode(origin) + threshold) / 2.;
        targetMiddleRight = (grid.getNode(targetEnd) + frac * grid.getNode(origin + 1) + threshold) / 2.;
        fractionLeft = (grid.getNode(targetBegin + 1) - (frac * grid.getNode(origin) + threshold)) / grid.duCell(origin);
        fractionRight = (frac * grid.getNode(origin + 1) + threshold - grid.getNode(targetEnd)) / grid.duCell(origin);
    } else
    {
        targetMiddleLeft = (grid.getNode(targetBegin + 1) + frac * grid.getNode(origin) - threshold) / 2.;
        targetMiddleRight = (grid.getNode(targetEnd) + frac * grid.getNode(origin + 1) - threshold) / 2.;
        fractionLeft = (grid.getNode(targetBegin + 1) - (frac * grid.getNode(origin) - threshold)) / grid.duCell(origin);
        fractionRight = (frac * grid.getNode(origin + 1) - threshold - grid.getNode(targetEnd)) / grid.duCell(origin);
    }

    double alphaMiddle = (grid.getNode(targetEnd) - grid.getNode(targetBegin + 1)) / grid.duCell(origin);
    double targetMiddleCenter = (grid.getNode(targetBegin + 1) + grid.getNode(targetEnd)) / 2.;

    auto alphaMin = alphaDistribution(targetMiddleLeft, grid.getCell(targetBegin), targetMiddleCenter, fractionLeft);
    alphaMiddle += alphaMin[1];
    auto alphaPlus = alphaDistribution(targetMiddleRight, targetMiddleCenter, grid.getCell(targetEnd), fractionRight);
    alphaMiddle += alphaPlus[0];

    std::vector<std::tuple<int, double>> alpha;
    alpha.push_back(std::make_tuple(targetBegin, alphaMin[0]));
    for (Grid::Index i = targetBegin + 1; i < targetEnd; i++)
    {
        alpha.push_back(std::make_tuple(i, grid.duCell(i)/(grid.getNode(targetEnd) - grid.getNode(targetBegin + 1)) * alphaMiddle));
    }
    alpha.push_back(std::make_tuple(targetEnd, alphaPlus[1]));

    return alpha;
}


std::vector<std::tuple<int, double>> getOperatorDistribution(const Grid& grid, double threshold, double source, int sourceidx, 
                                                             bool reverse = 0, double frac = 1.0)
{
    double targetCell;
    int targetBegin;
    int targetEnd;
    if (reverse)
    {
        targetCell = frac * source + threshold;
        targetBegin = std::upper_bound(grid.getNodes().begin(), grid.getNodes().end(), frac * grid.getNode(sourceidx) + threshold) - 
            grid.getNodes().begin() - 1;
        targetEnd = std::upper_bound(grid.getNodes().begin() + targetBegin, grid.getNodes().end(), frac * grid.getNode(sourceidx + 1) + threshold) -
            grid.getNodes().begin() - 1;
    } else
    {
        targetCell = frac * source - threshold;
        targetBegin = std::upper_bound(grid.getNodes().begin(), grid.getNodes().end(), frac * grid.getNode(sourceidx) - threshold) - 
            grid.getNodes().begin() - 1;
        targetEnd = std::upper_bound(grid.getNodes().begin() + targetBegin, grid.getNodes().end(), frac * grid.getNode(sourceidx + 1) - threshold) -
            grid.getNodes().begin() - 1;
    }

    std::vector<std::tuple<int, double>> alpha;
    if (targetBegin == targetEnd)
    {
        alpha = distributeOneCell(grid, targetCell, targetBegin);
    } else if (targetEnd - targetBegin == 1)
    {
        alpha = distributeTwoCells(grid,targetCell, targetBegin, targetEnd);
    } else
    {
        alpha = distributeNCells(grid, targetCell, targetBegin, targetEnd, sourceidx, threshold, reverse);
    }

    return alpha;
}

std::vector<std::tuple<int, double>> distributeNCellsIonization(const Grid& grid, double targetCell, int targetBegin, int targetEnd, 
    Grid::Index origin)
{
    double targetMiddleRight;
    double fractionRight;
    targetMiddleRight = (grid.getNode(targetEnd) + grid.getNode(origin + 1) - grid.getNode(origin)) / 2.;
    fractionRight = (grid.getNode(origin + 1) - grid.getNode(origin) - grid.getNode(targetEnd)) / grid.duCell(origin);

    double alphaMiddle = (grid.getNode(targetEnd) - grid.getNode(targetBegin)) / grid.duCell(origin);
    double targetMiddleCenter = (grid.getNode(targetBegin) + grid.getNode(targetEnd)) / 2.;

    auto alphaPlus = alphaDistribution(targetMiddleRight, targetMiddleCenter, grid.getCell(targetEnd), fractionRight);
    alphaMiddle += alphaPlus[0];

    std::vector<std::tuple<int, double>> alpha;
    for (Grid::Index i = targetBegin; i < targetEnd; i++)
    {
        alpha.push_back(std::make_tuple(i, grid.duCell(i)/(grid.getNode(targetEnd) - grid.getNode(targetBegin)) * alphaMiddle));
    }
    alpha.push_back(std::make_tuple(targetEnd, alphaPlus[1]));

    return alpha;
}

std::vector<std::tuple<int, double>> oneTakesAllDistribution(const Grid& grid, int sourceidx)
{
    double targetCell;
    int targetBegin;
    int targetEnd;
    targetCell = grid.getCell(sourceidx) - grid.getNode(sourceidx);
    targetBegin = 0;
    targetEnd = std::upper_bound(grid.getNodes().begin(), grid.getNodes().end(), grid.getNode(sourceidx + 1) - grid.getNode(sourceidx)) -
        grid.getNodes().begin() - 1;

    if (targetEnd == 0)
    {
        std::vector<std::tuple<int, double>> alpha;
        alpha.push_back(std::make_tuple(0, 1.));
        return alpha;
    }

    std::vector<std::tuple<int, double>> alpha;
    if (targetBegin == targetEnd)
    {
        alpha = distributeOneCell(grid, targetCell, targetBegin);
    } else if (targetEnd - targetBegin == 1)
    {
        alpha = distributeTwoCells(grid,targetCell, targetBegin, targetEnd);
    } else
    {
        alpha = distributeNCellsIonization(grid, targetCell, targetBegin, targetEnd, sourceidx);
    }

    return alpha;
}

InelasticOperator::InelasticOperator(const Grid& grid)
: hasSuperelastics(false)
{
    inelasticMatrix.setZero(grid.nCells(), grid.nCells());
}

void InelasticOperator::evaluateInelasticOperators(const Grid& grid, const EedfMixture& mixture)
{
    const Grid::Index cellNumber = grid.nCells();
    inelasticMatrix.setZero();
    /** \todo the following line seems to be missing. See the notes in the
     *  class declaration of this member
     */
    // hasSuperelastics = false;

    for (const auto &cd : mixture.collision_data().data_per_gas())
    {
        for (auto vecIndex : { CollisionType::excitation, CollisionType::vibrational, CollisionType::rotational})
        {
            for (const auto &collision : cd.collisions(vecIndex) )
            {
                const double threshold = collision->crossSection->threshold();

                if (threshold < grid.getNode(0) || threshold > grid.getNode(cellNumber))
                    continue;

                const double targetDensity = collision->getTarget()->delta();

                if (targetDensity != 0)
                {
                    Vector cellCrossSection(cellNumber);
                    Grid::Index numThreshold;

                    for (Grid::Index i = 0; i < cellNumber; ++i)
                        cellCrossSection[i] = 0.5 * ((*collision->crossSection)[i] + (*collision->crossSection)[i + 1]);

                    if (grid.isUniform())
                    {
                        numThreshold = static_cast<Grid::Index>(std::floor(threshold / grid.du()));
                    } else
                    {
                       numThreshold = std::upper_bound(grid.getNodes().begin(),grid.getNodes().end(), threshold) - grid.getNodes().begin() - 1;
                    }

                    if (grid.isUniform())
                    {
                        for (Grid::Index k = 0; k < cellNumber; ++k)
                        {
                            if (k < cellNumber - numThreshold)
                                inelasticMatrix(k, k + numThreshold) +=
                                    targetDensity * grid.getCells()[k + numThreshold] * cellCrossSection[k + numThreshold];

                            /** \todo Clarify. See the comments on the (conserving) attachment operator.
                             */
                            inelasticMatrix(k, k) -= targetDensity * grid.getCell(k) * cellCrossSection[k];
                        }
                    } else
                    {
                        for (Grid::Index k = 0; k < cellNumber; ++k)
                        {
                            if (grid.getNode(k+1) + threshold < grid.getCell(grid.nCells() - 1))\
                            {
                                std::vector<std::tuple<int, double>> alpha = getOperatorDistribution(grid, threshold,
                                     grid.getCell(k), k, true);

                                for (int i = 0; i < int(alpha.size()); i++)
                                {
                                     inelasticMatrix(k, std::get<0>(alpha[i])) += std::get<1>(alpha[i]) *
                                        targetDensity * grid.getCells()[std::get<0>(alpha[i])] * cellCrossSection[std::get<0>(alpha[i])];
                                }
                            }
                            inelasticMatrix(k, k) -= targetDensity * grid.getCell(k) * cellCrossSection[k];
                        }
                    }

                    if (collision->isReverse())
                    {
                        const double swRatio = collision->getTarget()->statisticalWeight /
                                               collision->m_rhsHeavyStates[0]->statisticalWeight;
                        const double productDensity = collision->m_rhsHeavyStates[0]->delta();

                        if (productDensity == 0)
                            continue;

                        /** \todo see the comments about superElasticThresholds in the ElectronKinetics.h.
                        if (numThreshold > 1)
                            superElasticThresholds.emplace_back(numThreshold);
                        */

                        if (numThreshold != 1)
                            hasSuperelastics = true;

                        if (grid.isUniform())
                        {
                            for (Grid::Index k = 0; k < cellNumber; ++k)
                            {
                                if (k >= numThreshold)
                                    inelasticMatrix(k, k - numThreshold) +=
                                        swRatio * productDensity * grid.getCell(k) * cellCrossSection[k];

                                if (k < cellNumber - numThreshold)
                                    inelasticMatrix(k, k) -= swRatio * productDensity * grid.getCell(k + numThreshold) *
                                                            cellCrossSection[k + numThreshold];
                            }
                        } else
                        {
                            for (Grid::Index k = 0; k < cellNumber; ++k)
                            {
                                if (grid.getNode(k) >= threshold)
                                {
                                    std::vector<std::tuple<int, double>> alpha = getOperatorDistribution(grid, threshold,
                                         grid.getCell(k), k);
                                    for (int i = 0; i < int(alpha.size()); i++)
                                    {
                                        inelasticMatrix(k, std::get<0>(alpha[i])) += std::get<1>(alpha[i]) *
                                            swRatio * productDensity * grid.getCells()[k] * cellCrossSection[k];
                                    }
                                }

                                if (grid.getNode(k + 1) + threshold < grid.getCell(cellNumber - 1))
                                {

                                    std::vector<std::tuple<int, double>> alpha = getOperatorDistribution(grid, threshold,
                                         grid.getCell(k), k, true);

                                    for (int i = 0; i < int(alpha.size()); i++)
                                    {
                                        inelasticMatrix(k, k) -= std::get<1>(alpha[i]) *
                                            swRatio * productDensity * grid.getCells()[std::get<0>(alpha[i])] * cellCrossSection[std::get<0>(alpha[i])];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

ElectronElectronOperator::ElectronElectronOperator(const Grid& grid)
{
    initialize(grid);
}

void ElectronElectronOperator::initialize(const Grid& grid)
{
    m_g_ee=0.0;
    updateABMatrices(grid);
    m_A.setZero(grid.nCells());
    m_B.setZero(grid.nCells());
}

void ElectronElectronOperator::clear()
{
    m_g_ee = 0.;
    m_A.setZero();
    m_B.setZero();
}

void ElectronElectronOperator::updateABMatrices(const Grid& grid)
{
    if (!grid.isUniform())
    {
        throw std::runtime_error("ElectronElectronOperator does not support nonuniform grids.");
    }
    
    /* Note that not all elements get a value below, the others must be set to
     * zero. We achieve that by clearing the entire matrix first.
     */
    m_a.setZero(grid.nCells(), grid.nCells());

    /* 1. Calculate a_ij/g_ee (equation 39f from the 2_2_0 manual).
     * This is the part right of the curly brace, but *with* the minus
     * sign that appears at the start of the expression. NOTE that in
     * the MATLAB code also a factor du is applied, which is removed
     * again when A and B are calculated from a_ij. In the C++ code
     * the matrix is really a_ij/g.
     */

    const Vector cellsThreeOverTwo = grid.getCells().cwiseProduct(grid.getCells().cwiseSqrt());
    const Vector energyArray = -(1./2.) * grid.getCells().cwiseSqrt() + (2./3./grid.du()) * cellsThreeOverTwo;
    for (Grid::Index i = 0; i < grid.nCells() - 1; ++i)
    {
        for (Grid::Index j = 1; j <= i; ++j)
        {
            m_a(i,j) = energyArray[j];
        }
        const double tmp = 2.0/3.0*std::pow(grid.getNodes()(i+1),1.5)/grid.du();
        for (Grid::Index j = i+1; j < grid.nCells(); ++j)
        {
            m_a(i,j) = tmp;
        }
    }

#define LOKIB_EE_APPLY_DB_FIX 1
#if LOKIB_EE_APPLY_DB_FIX == 1
    const Matrix tmp(m_a);
    for (Grid::Index i = 0; i < grid.nCells() - 1; ++i)
    {
        for (Grid::Index j = 1; j < grid.nCells(); ++j)
        {
            m_a(i,j) = std::sqrt(tmp(i,j) * tmp(j-1,i+1));
        }
    }
#else
    for (Grid::Index i = 0; i < grid.nCells() - 1; ++i)
    {
        for (Grid::Index j = 1; j < grid.nCells(); ++j)
        {
            m_a(i,j) = std::sqrt(m_a(i,j) * m_a(j-1,i+1));
        }
    }
#endif // LOKIB_EE_APPLY_DB_FIX
}

void ElectronElectronOperator::update_g_ee_AB(const Grid& grid, const Vector& eedf, double ne, double n0)
{
    const double e = Constant::electronCharge;
    const double e0 = Constant::vacuumPermittivity;
    const Vector cellsThreeOverTwo = grid.getCells().cwiseProduct(grid.getCells().cwiseSqrt());

    const double meanEnergy = grid.du() * cellsThreeOverTwo.dot(eedf);
    const double Te = 2. / 3. * meanEnergy;
    const double logC = std::log(12 * Constant::pi * std::pow(e0 * Te / e, 1.5) / std::sqrt(ne));
    // g_ee = (3*alpha_eV/gamma)*(ne/N)
    m_g_ee = (ne / n0) * (e * e / (8 * Constant::pi * e0 * e0)) * logC;
    m_A = m_g_ee * (m_a * eedf);
    m_B = m_g_ee * (m_a.transpose() * eedf);
}

void ElectronElectronOperator::evaluatePower(const Grid& grid, const Vector& eedf, double& power) const
{
    /*  One du() comes from the integration, together with eedf).
     *  The other is the energy gain of an electron moving up one cell,
     *  I (JvD) believe. Make sure this is explained well in the docs.
     */
    power = (-SI::gamma * grid.du() * grid.du()) * (m_A - m_B).dot(eedf);
    //std::cout << "EE POWER: " << power << std::endl;
}

void ElectronElectronOperator::discretizeTerm(Matrix& M, const Grid& grid) const
{
    if (!grid.isUniform())
    {
        throw std::runtime_error("ElectronElectronOperator does not support nonuniform grids.");
    }
    
    for (Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        M(k,k) += - (m_A[k] + m_B[k]);
        if (k > 0)
        {
            M(k,k-1) += m_A[k - 1];
        }
        if (k < grid.nCells() - 1)
        {
            M(k,k+1) += m_B[k + 1];
        }
    }
}

IonizationOperator::IonizationOperator(IonizationOperatorType type)
: ionizationOperatorType(type), includeNonConservativeIonization(false)
{
}

void IonizationOperator::evaluateIonizationOperator(const Grid& grid, const EedfMixture& mixture)
{
    bool hasValidCollisions = false;

    ionConservativeMatrix.setZero();

    if (ionizationOperatorType != IonizationOperatorType::conservative)
        ionizationMatrix.setZero();
    /** \todo the following line seems to be missing. See the notes in the
     *  class declaration of this member
     */
    //includeNonConservativeIonization = false;

    for (const auto &cd : mixture.collision_data().data_per_gas())
    {
        for (const auto &collision : cd.collisions(CollisionType::ionization))
        {
            const double threshold = collision->crossSection->threshold();

            if (threshold > grid.getNode(grid.nCells()))
                continue;

            hasValidCollisions = true;
            const double delta = collision->getTarget()->delta();

            Vector cellCrossSection(grid.nCells());
            collision->crossSection->interpolate(grid.getCells(), cellCrossSection);
            int numThreshold1 = static_cast<uint32_t>(std::upper_bound(grid.getCells().begin(),grid.getCells().end(), threshold) - grid.getCells().begin() - 1);
            int numThreshold2 = numThreshold1 + 1;
            int numThreshold = std::abs(threshold - grid.getCell(numThreshold1)) < std::abs(threshold - grid.getCell(numThreshold2)) ? numThreshold1 : numThreshold2;
           
            switch (ionizationOperatorType)
            {
            case IonizationOperatorType::conservative:
                break;

            case IonizationOperatorType::oneTakesAll:
                if (grid.isUniform())
                {
                    for (Grid::Index k = 0; k < grid.nCells(); ++k)
                    {
                        // Manual_2_2_0 eq. 15. TODO: Note that newborns are
                        // inserted in the first cell, so at du/2, not at zero.
                        if (k < grid.nCells() - numThreshold)
                            ionizationMatrix(k, k + numThreshold) +=
                                delta * grid.getCell(k + numThreshold) * cellCrossSection[k + numThreshold];

                        const double term = delta * grid.getCell(k) * cellCrossSection(k);

                        ionizationMatrix(k, k) -= term;
                        ionizationMatrix(0, k) += term;
                    }
                } else
                {
                    for (Grid::Index k = 0; k < grid.nCells(); ++k)
                    {
                        // Manual_2_2_0 eq. 15. TODO: Note that newborns are
                        // inserted in the first cell, so at du/2, not at zero.
                        if (grid.getNode(k+1) + threshold < grid.getCell(grid.nCells() - 1))
                        {
                            std::vector<std::tuple<int, double>> alpha = getOperatorDistribution(grid, threshold,
                                    grid.getCell(k), k, true);

                            for (int i = 0; i < int(alpha.size()); i++)
                            {
                                ionizationMatrix(k, std::get<0>(alpha[i])) += std::get<1>(alpha[i]) *
                                    delta * grid.getCell(std::get<0>(alpha[i])) * cellCrossSection[std::get<0>(alpha[i])];
                            }
                        }

                        const double term = delta * grid.getCell(k) * cellCrossSection(k);

                        ionizationMatrix(k, k) -= term;

                        std::vector<std::tuple<int, double>> alpha = oneTakesAllDistribution(grid, k);

                        for (int i = 0; i < int(alpha.size()); i++)
                        {
                            ionizationMatrix(std::get<0>(alpha[i]), k) += std::get<1>(alpha[i]) * term;
                        }
                    }
                }
                
                break;

            case IonizationOperatorType::equalSharing:
                if (grid.isUniform())
                {
                    for (Grid::Index k = 0; k < grid.nCells(); ++k)
                    {
                        // Manual_2_2_0 eq. 16
                        ionizationMatrix(k, k) -= delta * grid.getCell(k) * cellCrossSection[k];

                        if (k < (grid.nCells() - numThreshold) / 2)
                        {
                            const Grid::Index i = 2 * (k + 1) + numThreshold - 1;

                            ionizationMatrix(k, i) += 4 * delta * grid.getCell(i) * cellCrossSection(i);
                        }
                    }
                } else
                {
                    for (Grid::Index k = 0; k < grid.nCells(); ++k)
                    {
                        ionizationMatrix(k, k) -= delta * grid.getCell(k) * cellCrossSection[k];

                        if (grid.getNode(k+1) < (grid.getCell(grid.nCells() - 1) - threshold) / 2)
                        {
                            std::vector<std::tuple<int, double>> alpha = getOperatorDistribution(grid, threshold,
                                    grid.getCell(k), k, true, 2.0);
                            for (int i = 0; i < int(alpha.size()); i++)
                            {
                                ionizationMatrix(k, std::get<0>(alpha[i])) += 4 * std::get<1>(alpha[i]) *
                                    delta * grid.getCell(std::get<0>(alpha[i])) * cellCrossSection[std::get<0>(alpha[i])];
                            }
                        }
                    }
                }
                
                break;
            case IonizationOperatorType::sdcs:

                // This case is in eq. 13b of Manual_2_2_0
                /// \todo Provide the missing details. Part of the discussion is around eq. 63.
                double W = cd.OPBParameter();

                if (W < 0)
                    W = threshold;

                for (Grid::Index k = 0; k < grid.nCells(); ++k)
                {
                    const Grid::Index end = std::min(2 * (k + 1) + numThreshold, grid.nCells());

                    if (k > numThreshold)
                    {
                        const Grid::Index half = (k + 1 - numThreshold) / 2;
                        const double numerator = 1 / std::atan((grid.getCell(k) - threshold) / (2 * W));
                        double sum = 0.;

                        for (Grid::Index i = 0; i < half; ++i)
                            sum += numerator / (W + grid.getCell(i) * grid.getCell(i) / W);

                        if (grid.isUniform())
                        {
                            ionizationMatrix(k, k) -= delta * grid.du() * grid.getCell(k) * cellCrossSection[k] * sum;
                        } else
                        {
                            ionizationMatrix(k, k) -= delta * grid.duCell(k) * grid.getCell(k) * cellCrossSection[k] * sum;
                        }
                    }

                    /** \todo If k + numThreshold + 1 >= grid.nCells(), the term is ignored.
                     *        Document (in the document, not necessarily here) what are the
                     *        consequences of that.
                     */
                    if (k + numThreshold + 1 < grid.nCells())
                    {
                        if (grid.isUniform())
                        {
                            for (Grid::Index i = k + numThreshold + 1; i < end; ++i)
                            {
                                ionizationMatrix(k, i) += delta * grid.du() * grid.getCell(i) * cellCrossSection[i] /
                                                        (std::atan((grid.getCell(i) - threshold) / (2 * W)) *
                                                        (W + std::pow(grid.getCell(i - k - numThreshold - 1), 2) / W));
                            /* NOTE: The expression looks a bit different from eq. 63 of the manual (2.2.0)
                            * because the w from the denominator of the first multiplicand is combined
                            * with the factor 1+(u/w)^beta, and beta=2 is hardcoded, resulting in the
                            * factor 1/(w+u^2/w).
                            */
                            }
                        } else 
                        {
                            for (Grid::Index i = k + numThreshold + 1; i < end; ++i)
                            {
                                ionizationMatrix(k, i) += delta * grid.duCell(i) * grid.getCell(i) * cellCrossSection[i] /
                                                        (std::atan((grid.getCell(i) - threshold) / (2 * W)) *
                                                        (W + std::pow(grid.getCell(i - k - numThreshold - 1), 2) / W));
                            }
                        }
                    }

                    /** \todo The following comment needs to be sorted out (possible index errors).
                     *  \todo Document what is done here (algorithm) and how it is implemented.
                     */

                    /** \todo This last section might need some adjustments because of indexing
                     *  differences between Matlab and C++ (since indexes are multiplied here).
                     */

                    for (Grid::Index i = 2 * (k + 1) + numThreshold - 1; i < grid.nCells(); ++i)
                    {
                        if (grid.isUniform())
                        {
                            ionizationMatrix(k, i) += delta * grid.du() * grid.getCell(i) * cellCrossSection[i] /
                                                    (std::atan((grid.getCell(i) - threshold) / (2 * W)) *
                                                    (W + std::pow(grid.getCell(k), 2) / W));
                            /* See the note above about the optical differences with equation 63
                            * of the manual (2.2.0).
                            */
                        } else
                        {
                            ionizationMatrix(k, i) += delta * grid.duCell(i) * grid.getCell(i) * cellCrossSection[i] /
                                                    (std::atan((grid.getCell(i) - threshold) / (2 * W)) *
                                                    (W + std::pow(grid.getCell(k), 2) / W));
                        }
                    }
                }
                break;
            }

            // Evaluation of the conservative ionization operator

            if (numThreshold == 0)
                continue;

            if (grid.isUniform())
            {
                for (Grid::Index k = 0; k < grid.nCells(); ++k)
                {
                    if (k < grid.nCells() - numThreshold)
                        ionConservativeMatrix(k, k + numThreshold) +=
                            delta * grid.getCell(k + numThreshold) * cellCrossSection[k + numThreshold];

                    ionConservativeMatrix(k, k) -= delta * grid.getCell(k) * cellCrossSection[k];
                }
            } else
            {
                for (Grid::Index k = 0; k < grid.nCells(); ++k)
                {
                    if (grid.getNode(k+1) + threshold < grid.getCell(grid.nCells() - 1))
                    {
                        std::vector<std::tuple<int, double>> alpha = getOperatorDistribution(grid, threshold,
                                grid.getCell(k), k, true);

                        for (int i = 0; i < int(alpha.size()); i++)
                        {
                                ionConservativeMatrix(k, std::get<0>(alpha[i])) += std::get<1>(alpha[i]) *
                                delta * grid.getCell(std::get<0>(alpha[i])) * cellCrossSection[std::get<0>(alpha[i])];
                        }
                    }
                    ionConservativeMatrix(k, k) -= delta * grid.getCell(k) * cellCrossSection[k];
                }
            }
            

            /// \todo Should this line not appear before the 'if (numThreshold == 0) continue;' test?
            if (ionizationOperatorType != IonizationOperatorType::conservative && hasValidCollisions)
                includeNonConservativeIonization = true;
        }
    }
}

AttachmentOperator::AttachmentOperator()
 : includeNonConservativeAttachment(false)
{
}

void AttachmentOperator::evaluateAttachmentOperator(const Grid& grid, const EedfMixture& mixture)
{

    attachmentMatrix.setZero();
    attachmentConservativeMatrix.setZero();

    const Grid::Index cellNumber = grid.nCells();

    for (const auto &cd : mixture.collision_data().data_per_gas())
    {
        for (const auto &collision : cd.collisions(CollisionType::attachment))
        {
            const double threshold = collision->crossSection->threshold();

            if (threshold > grid.getNode(cellNumber))
                continue;

            /* This should definitely not be in this (double) loop. Is this a constructor task?
             * Can this just be replaced with 'cd.collisions(CollisionType::attachment).size()'
             * in (other) places where this is now used?
             * Answer: no, this depends on uMax(), which may change for a smart grid. But it
             * could be done as an action when that changes.
             */
            /** \bug This should be reset to false at the beginning of this function
             *       because the results may change when uMax is changed.
             */
            includeNonConservativeAttachment = true;

            /** \todo Eliminate the cellCrossSection vector? This can be calculated on the fly
             *        in the two places where it is needed (one if the merger below can be done).
             */
            Vector cellCrossSection(cellNumber);

            const double targetDensity = collision->getTarget()->delta();

            /// \todo Merge with the subsequent k-loop.
            for (Grid::Index i = 0; i < cellNumber; ++i)
                cellCrossSection[i] = 0.5 * ((*collision->crossSection)[i] + (*collision->crossSection)[i + 1]);

            // This is eq. 13c of \cite Manual_2_2_0
            for (Grid::Index k = 0; k < cellNumber; ++k)
                attachmentMatrix.coeffRef(k, k) -= targetDensity * grid.getCell(k) * cellCrossSection[k];

            int numThreshold;
            if (grid.isUniform())
            {
                numThreshold = static_cast<Grid::Index>(std::floor(threshold / grid.du()));
            } else
            {
                Grid::Index findIndex = (std::upper_bound(grid.getNodes().begin(),grid.getNodes().end(), threshold) - grid.getNodes().begin());
                numThreshold = static_cast<Grid::Index>(findIndex) - 1;
            }
            /* This is an optimization: do not add a source and a sink that cancel out.
             * When this happens, the (minor) energy gain is ignored.
             */
            if (numThreshold == 0)
                continue;

            /** Can we also merge with this loop?
             *  It does not seem problematic to have an 'if (numThreshold)' in the loop,
             *  given the amount of work done.
             */
            for (Grid::Index k = 0; k < cellNumber; ++k)
            {
                /** \todo The manual \cite Manual_2_2_0 states that attachment is always
                 *  non-conserving (5 lines below eq. 16). What is attachmentConservativeMatrix,
                 *  then?
                 */
                /** \todo Explain the structure of this term. Explain that this is conserving,
                 *  how we can see that (integral source must be zero).
                 */
                /** \todo Is this really a conserving term? It appears that this loses electrons
                 *  if the argument of the following if-statement is false (that is: for
                 *  electrons with at a distance smaller than the attachment energy from uMax.
                 *  Of course this is in practice only a minor problem. It can be fixed by
                 *  also having the sink inside the if-statement, so we simply discard the
                 *  processes that produce electrons with an energy that cannot be represented.
                 */
                if (k < cellNumber - numThreshold)
                    attachmentConservativeMatrix(k, k + numThreshold) +=
                        targetDensity * grid.getCell(k + numThreshold) * cellCrossSection[k + numThreshold];

                attachmentConservativeMatrix(k, k) -= targetDensity * grid.getCell(k) * cellCrossSection[k];
            }
        }
    }
}

} // namespace loki
