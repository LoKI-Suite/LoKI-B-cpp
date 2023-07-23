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

InelasticOperator::InelasticOperator(const Grid& grid)
: hasSuperelastics(false)
{
    inelasticMatrix.setZero(grid.nCells(), grid.nCells());
}

void InelasticOperator::evaluateInelasticOperators(const Grid& grid, const EedfMixture& mixture)
{
    const Grid::Index cellNumber = grid.nCells();

    inelasticMatrix.setZero();

    for (const auto &cd : mixture.collision_data().data_per_gas())
    {
        for (auto vecIndex : { CollisionType::excitation, CollisionType::vibrational, CollisionType::rotational})
        {
            for (const auto &collision : cd.collisions(vecIndex) )
            {
                const double threshold = collision->crossSection->threshold();

                if (threshold < grid.du() || threshold > grid.getNodes()[grid.nCells()])
                    continue;

                const double targetDensity = collision->getTarget()->delta();

                if (targetDensity != 0)
                {
                    const auto numThreshold = static_cast<Grid::Index>(std::floor(threshold / grid.du()));

                    Vector cellCrossSection(cellNumber);

                    for (Grid::Index i = 0; i < cellNumber; ++i)
                        cellCrossSection[i] = 0.5 * ((*collision->crossSection)[i] + (*collision->crossSection)[i + 1]);

                    for (Grid::Index k = 0; k < cellNumber; ++k)
                    {
                        if (k < cellNumber - numThreshold)
                            inelasticMatrix(k, k + numThreshold) +=
                                targetDensity * grid.getCells()[k + numThreshold] * cellCrossSection[k + numThreshold];
                        /** \todo Clarify. See the comments on the (conserving) attachment operator.
                         */
                        inelasticMatrix(k, k) -= targetDensity * grid.getCells()[k] * cellCrossSection[k];
                    }

                    if (collision->isReverse())
                    {
                        const double swRatio = collision->getTarget()->statisticalWeight /
                                               collision->m_rhsHeavyStates[0]->statisticalWeight;
                        const double productDensity = collision->m_rhsHeavyStates[0]->delta();

                        if (productDensity == 0)
                            continue;

                        /** \todo see the comments about superElasticThresholds in the header file.
                        if (numThreshold > 1)
                            superElasticThresholds.emplace_back(numThreshold);
                        */

                        if (numThreshold != 1)
                            hasSuperelastics = true;

                        for (Grid::Index k = 0; k < cellNumber; ++k)
                        {
                            if (k >= numThreshold)
                                inelasticMatrix(k, k - numThreshold) +=
                                    swRatio * productDensity * grid.getCells()[k] * cellCrossSection[k];

                            if (k < cellNumber - numThreshold)
                                inelasticMatrix(k, k) -= swRatio * productDensity * grid.getCells()[k + numThreshold] *
                                                         cellCrossSection[k + numThreshold];
                        }
                    }
                }
            }
        }
    }
}

ElectronElectronOperator::ElectronElectronOperator()
{
    g_ee=0.0;
}

void ElectronElectronOperator::initialize(const Grid& grid)
{
    A.setZero(grid.nCells());
    B.setZero(grid.nCells());
}

void ElectronElectronOperator::updateABMatrices(const Grid& grid)
{
    BAee.setZero(grid.nCells(), grid.nCells());

    /* 1. Calculate a_nm/g_ee (equation 39f from the 2_2_0 manual).
     * This is the part right of the curly brace, but with the minus sign
     * that appears at the start of the expression added.
     */

    const Vector cellsThreeOverTwo = grid.getCells().cwiseProduct(grid.getCells().cwiseSqrt());
    const Vector energyArray = -(grid.du() / 2.) * grid.getCells().cwiseSqrt() + (2. / 3.) * cellsThreeOverTwo;

    /** \todo It appears that we are setting up a^T here, instead of a.
     *  You could also say that we are setting up b here, since b = a^T
     *  is assumed later on. But just setting up a first seems to bring
     *  this closer to the documentation.
     *  \todo change the name. Just call this 'a', instead of 'BAee',
     *  now that this has become a member of ElectronElectronOperator.
     */
    for (Grid::Index j = 0; j < grid.nCells() - 1; ++j)
    {

        for (Grid::Index i = 1; i <= j; ++i)
            BAee(i, j) = energyArray[i];

        const double value = 2. / 3. * std::pow(grid.getNode(j + 1), 1.5);
        
        for (Grid::Index i = j + 1; i < grid.nCells(); ++i)
            BAee(i, j) = value;
    }

    // detailed balance condition

    for (Grid::Index j = 0; j < grid.nCells() - 1; ++j)
    {
        for (Grid::Index i = 1; i < grid.nCells(); ++i)
        {
            BAee(i, j) = std::sqrt(BAee(i, j) * BAee(j + 1, i - 1));
        }
    }
}

void ElectronElectronOperator::update_g_ee(const Grid& grid, const Vector& eedf, double ne, double n0)
{
    const double e = Constant::electronCharge;
    const double e0 = Constant::vacuumPermittivity;
    const Vector cellsThreeOverTwo = grid.getCells().cwiseProduct(grid.getCells().cwiseSqrt());

    double meanEnergy = grid.du() * cellsThreeOverTwo.dot(eedf);
    double Te = 2. / 3. * meanEnergy;
    double logC = std::log(12 * Constant::pi * std::pow(e0 * Te / e, 1.5) / std::sqrt(ne));
    g_ee = (ne / n0) * (e * e / (8 * Constant::pi * e0 * e0)) * logC;
}

void ElectronElectronOperator::updateAB(const Grid& grid, const Vector& eedf)
{
        /** \todo Should .transpose() not be applied to the 2nd expression instead of the first?
         *  Also see the note about a^T vs. a above. That explains the lines below. Is there a
         *  reason to do things this way?
         */
        A = (g_ee / grid.du()) * (BAee.transpose() * eedf);
        B = (g_ee / grid.du()) * (BAee * eedf);
}

void ElectronElectronOperator::evaluatePower(const Grid& grid, const Vector& eedf, double& power) const
{
        /** \todo One du() comes from the integration, together with eedf).
         *  It feels odd that the second is not part of the definitions of A and B:
         *  check why this is not the case.
         */
        power = (-SI::gamma * grid.du() * grid.du()) * (A - B).dot(eedf);
        //std::cout << "EE POWER: " << power << std::endl;
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
            const auto numThreshold = static_cast<Grid::Index>(std::floor(threshold / grid.du()));

            Vector cellCrossSection(grid.nCells());

            for (Grid::Index i = 0; i < grid.nCells(); ++i)
                cellCrossSection[i] = 0.5 * ((*collision->crossSection)[i] + (*collision->crossSection)[i + 1]);

            switch (ionizationOperatorType)
            {
            case IonizationOperatorType::conservative:
                break;

            case IonizationOperatorType::oneTakesAll:
                for (Grid::Index k = 0; k < grid.nCells(); ++k)
                {
                    if (k < grid.nCells() - numThreshold)
                        ionizationMatrix(k, k + numThreshold) +=
                            delta * grid.getCell(k + numThreshold) * cellCrossSection[k + numThreshold];

                    const double term = delta * grid.getCell(k) * cellCrossSection(k);

                    ionizationMatrix(k, k) -= term;
                    ionizationMatrix(0, k) += term;
                }
                break;

            case IonizationOperatorType::equalSharing:
                for (Grid::Index k = 0; k < grid.nCells(); ++k)
                {
                    ionizationMatrix(k, k) -= delta * grid.getCell(k) * cellCrossSection[k];

                    if (k < (grid.nCells() - numThreshold) / 2)
                    {
                        const Grid::Index i = 2 * (k + 1) + numThreshold - 1;

                        ionizationMatrix(k, i) += 4 * delta * grid.getCell(i) * cellCrossSection(i);
                    }
                }
                break;
            case IonizationOperatorType::sdcs:
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

                        ionizationMatrix(k, k) -= delta * grid.du() * grid.getCell(k) * cellCrossSection[k] * sum;
                    }

                    /** \todo If k + numThreshold + 1 >= grid.nCells(), the term is ignored.
                     *        Document (in the document, not necessarily here) what are the
                     *        consequences of that.
                     */
                    if (k + numThreshold + 1 < grid.nCells())
                    {
                        for (Grid::Index i = k + numThreshold + 1; i < end; ++i)
                        {
                            ionizationMatrix(k, i) += delta * grid.du() * grid.getCell(i) * cellCrossSection[i] /
                                                      (std::atan((grid.getCell(i) - threshold) / (2 * W)) *
                                                       (W + std::pow(grid.getCell(i - k - numThreshold - 1), 2) / W));
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
                        ionizationMatrix(k, i) += delta * grid.du() * grid.getCell(i) * cellCrossSection[i] /
                                                  (std::atan((grid.getCell(i) - threshold) / (2 * W)) *
                                                   (W + std::pow(grid.getCell(k), 2) / W));
                    }
                }
                break;
            }

            // Evaluation of the conservative ionization operator

            if (numThreshold == 0)
                continue;

            for (Grid::Index k = 0; k < grid.nCells(); ++k)
            {
                if (k < grid.nCells() - numThreshold)
                    ionConservativeMatrix(k, k + numThreshold) +=
                        delta * grid.getCell(k + numThreshold) * cellCrossSection[k + numThreshold];

                ionConservativeMatrix(k, k) -= delta * grid.getCell(k) * cellCrossSection[k];
            }

            /// \todo Should this line not appear before the 'if (numThreshold == 0) continue;' test?
            if (ionizationOperatorType != IonizationOperatorType::conservative && hasValidCollisions)
                includeNonConservativeIonization = true;
        }
    }
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

            for (Grid::Index k = 0; k < cellNumber; ++k)
                attachmentMatrix.coeffRef(k, k) -= targetDensity * grid.getCell(k) * cellCrossSection[k];

            const auto numThreshold = static_cast<Grid::Index>(std::floor(threshold / grid.du()));
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
                /** \todo Explain the structore of this term. Explain that this is conserving,
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
