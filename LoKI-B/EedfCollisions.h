/** \file
 *
 *  Interfaces of LoKI-B's code for electron-impact collisions and
 *  containers (per gas and for a mixture) of such collisions.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2020 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
 *  M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
 *  L. L. Alves
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  \author Daan Boer and Jan van Dijk (C++ version)
 *  \date   21. May 2019
 */

#ifndef LOKI_CPP_EEDFCOLLISIONS_H
#define LOKI_CPP_EEDFCOLLISIONS_H

#include "LoKI-B/Collision.h"
#include "LoKI-B/CrossSection.h"
#include "LoKI-B/Gas.h"
#include "LoKI-B/GasMixture.h"
#include "LoKI-B/Enumeration.h"
#include "LoKI-B/MacroscopicQuantities.h"
#include "LoKI-B/Grid.h"
#include "LoKI-B/Power.h"

#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

namespace loki
{

class EedfCollision : public Collision
{
public:
    using EedfState = Gas::State;

    EedfCollision(CollisionType type, const StateVector &lhsStates, const CoeffVector &lhsCoeffs,
                  const StateVector &rhsStates, const CoeffVector &rhsCoeffs, bool isReverse);
    ~EedfCollision();
    const EedfState *getTarget() const;
    EedfState *getTarget();
    void superElastic(const Vector &energyData, Vector &result) const;
    PowerTerm evaluateConservativePower(const Vector &eedf) const;
    PowerTerm evaluateNonConservativePower(const Vector &eedf, const IonizationOperatorType ionizationOperatorType,
                                           const double OPBParameter) const;
    /// \todo non-constant because ineRateCoeff, supRateCoeff are modified
    RateCoefficient evaluateRateCoefficient(const Vector &eedf);
    std::string typeAsString() const;
    friend std::ostream &operator<<(std::ostream &os, const EedfCollision &collision);

    double ineRateCoeff() const { return m_ineRateCoeff; }
    double supRateCoeff() const { return m_supRateCoeff; }
private:
    // The raw cross section data and threshold is stored in
    // the CrossSection object
    const StateVector m_lhsHeavyStates;
    const CoeffVector m_lhsHeavyCoeffs;
public:
    /// \todo Make private
    StateVector m_rhsHeavyStates;
private:
    CoeffVector m_rhsHeavyCoeffs;
    double m_ineRateCoeff;
    double m_supRateCoeff;
public:
    std::unique_ptr<CrossSection> crossSection;

};

// TODO: Allow loading of effective populations from a file.
class EedfCollisionDataGas
{
public:
    using State = Gas::State;
    using CollisionVector = std::vector<std::unique_ptr<EedfCollision>>;
    using CollisionsType = std::vector<CollisionVector>;
    EedfCollisionDataGas(const GasProperties& gasProps, Gas& gas);
    EedfCollisionDataGas(const EedfCollisionDataGas&) = delete;
    EedfCollisionDataGas(EedfCollisionDataGas&&) = default;
    EedfCollisionDataGas& operator=(const EedfCollisionDataGas&) = delete;
    EedfCollisionDataGas& operator=(EedfCollisionDataGas&&) = delete;
    ~EedfCollisionDataGas() {}

    // We need to store the collisions per Gas since we need to calculate
    // the mass ratio when evaluating the total and elastic cross-sections.
    const CollisionVector& collisions(CollisionType type) const { return m_collisions[static_cast<uint8_t>(type)]; }
    const CollisionsType& collisions() const { return m_collisions; }
    const CollisionsType& collisionsExtra() const { return m_collisionsExtra; }
    void addCollision(EedfCollision *collision, bool isExtra);
    void checkElasticCollisions(State *electron, const Grid *energyGrid);
    bool isDummy() const;
    const GasPower &getPower() const;
    /** \todo Non-constant because m_power is changed. See if m_power must managed here.
     *  The name updatePower would make more clear that this modifies the object.
     */
    const GasPower& evaluatePower(const IonizationOperatorType ionType, const Vector &eedf);
    double OPBParameter() const { return m_OPBParameter; }
    const Gas& gas() const { return m_gas; }
private:
    /* the following three members are used (only) for Effective -> Elastic,
     * (together with the public checkElasticCollisions).
     */
    void setDefaultEffPop(State *ground);
    CrossSection *elasticCrossSectionFromEffective(const Grid *energyGrid);
    std::vector<State *> findStatesToUpdate();
    // two helpers for evaluating power terms
    PowerTerm evaluateConservativePower(const CollisionVector &collisionVector, const Vector &eedf) const;
    PowerTerm evaluateNonConservativePower(const CollisionVector &collisionVector, const IonizationOperatorType ionType,
                                       const Vector &eedf) const;
    Gas& m_gas;
    CollisionsType m_collisions;
    CollisionsType m_collisionsExtra;
    std::map<const State *, std::vector<EedfCollision *>> m_state_collisions;
    std::map<const State *, std::vector<EedfCollision *>> m_state_collisionsExtra;
    std::map<const State *, double> m_effectivePopulations;
    GasPower m_power;
    double m_OPBParameter;
};

/** \todo See if we can somehow separate settings (fixed) and mutable
 *  output data in this class (and its members).
 */
class EedfCollisionDataMixture
{
public:
    using State = Gas::State;
    using Collision = EedfCollision;
    using EedfCollisionDataGasArray = std::vector<EedfCollisionDataGas>;
    struct CollisionEntry
    {
        Collision* m_coll;
        const bool m_isExtra;
    };
    EedfCollisionDataMixture();
    /** Loads the collisions from the \a file that is provided as first argument.
     *  It also needs a pointer to the \a energyGrid and a boolean \a isExtra to
     *  indicate whether the collisions are extra, for correct initialization and
     *  storage of the collisions.
     *  \todo Explain isExtra.
     */
    void loadCollisionsClassic(const std::filesystem::path &file, const GasProperties& gasProps, GasMixture& composition, const Grid *energyGrid, bool isExtra);

    /** Loads the collisions from a json mixture section. Such section must contain two
     *  subsections: an object "states" that consists of key-value pairs that represent
     *  name of the species and an object that defines the species, and an array called
     *  "sets". Ech set element is an object that describes a group of processes. It
     *  must contain an array "processes"; each element of this array describes a processes.
     *  Furthermore, it needs a pointer to the energy grid and a boolean to
     *  indicate whether the collisions are extra, for correct initialization and storage of
     *  the collisions.
     *
     *  \todo The meta information in the set objects ("_id", "complete", "description",
     *  "contributor" etc. are ignored so far.
     */
    void loadCollisionsJSON(const json_type &mcnf, const GasProperties& gasProps, GasMixture& composition, const Grid *energyGrid, bool isExtra);

    const EedfCollisionDataGasArray& data_per_gas() const { return m_data_per_gas; }
    EedfCollisionDataGasArray& data_per_gas() { return m_data_per_gas; }
    std::vector<CollisionEntry> m_collisions;
    const EedfCollisionDataGas& get_data_for(const Gas& gas) const;
    EedfCollisionDataGas& get_data_for(const Gas& gas);
    const Vector& elasticCrossSection() const { return m_elasticCrossSection; }
    const Vector& totalCrossSection() const { return m_totalCrossSection; }
    /** \todo Add precise references to published sources for this function.
     *  The elastic cross section implemented here is related to the
     *  combination of equations 12 and equation 6b of \cite Tejero2019
     *  and the definitions just below 6d in that paper. It appears that
     *  the exact expressions implemented in this function are not
     *  explicitly stated in the text. Not a problem per se, but we should
     *  make sure to compensate for that by having detailed documentation here.
     */
    void evaluateTotalAndElasticCS(const Grid& grid);
    /** \bug The semantics of this member are not clear. Should this consider the 'extra'
     *  collisions as well? At present, they are (in loadCollisions the true flag is set
     *  when a collision of a type is created). At the same time, hasCollisions is used
     *  in ElectronKinetics.cpp to decide whether to add particular terms to the Boltzmann
     *  matrix. If, subsequently, the collsions(that_type) vector is used to populate those
     *  matrix contributions, zero-valued terms may be added (if only *extra* collisions
     *  of that type were loaded).
     *  \todo Should hasCollisions be recalculated when uMax() is changed?
     *  (smartGrid)? In that case process types may coma and go if they
     *  have threshold above uMax().
     */
    bool hasCollisions(CollisionType type) const { return m_hasCollisions[static_cast<uint8_t>(type)]; }
    /// \todo make private, provide const accessors?
    std::vector<RateCoefficient> m_rateCoefficients;
    std::vector<RateCoefficient> m_rateCoefficientsExtra;
    /** \bug This ADDS rate coefficients to the previously calulated ones in
     *       rateCoefficients and rateCoefficientsExtra. Every time a case is run,
     *       the list grows. This is obviously wrong, but before changing this, let's
     *       see what the MATLAB version of LoKI-B does exactly.
     *  \todo Should processes with a too-high threshold be skipped, as is done now?
     *        This may result in tables with different sizes when multiple cases are run.
     *        It may be better to just produce zero-values entries in such case.
     *  \todo Also the rate coefficient for the effective cross section is produced.
     *        Is that a meaningful value?
     */
    void evaluateRateCoefficients(const Grid& grid, const Vector &eedf);
private:
    /// \todo Document me more precisely: returns nullptr for a duplicate proccess
    EedfCollision* addCollision(CollisionType type,
        const Collision::StateVector &lhsStates, const Collision::CoeffVector &lhsCoeffs,
        const Collision::StateVector &rhsStates, const Collision::CoeffVector &rhsCoeffs,
        bool reverseAlso, bool isExtra);
    State *ensureState(const GasProperties& gasProps, GasMixture& composition, const StateEntry &entry);
    EedfCollisionDataGasArray m_data_per_gas;
    bool m_hasCollisions[static_cast<uint8_t>(CollisionType::size)]{false};
    Vector m_elasticCrossSection;
    Vector m_totalCrossSection;
};

} // namespace loki

#endif // LOKI_CPP_EEDFCOLLISIONS_H
