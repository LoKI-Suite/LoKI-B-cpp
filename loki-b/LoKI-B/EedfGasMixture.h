//
// Created by daan on 20-5-19.
//

#ifndef LOKI_CPP_EEDFGASMIXTURE_H
#define LOKI_CPP_EEDFGASMIXTURE_H

#include "LoKI-B/EedfCollision.h"
#include "LoKI-B/EedfGas.h"
#include "LoKI-B/GasMixture.h"
#include "LoKI-B/MacroscopicQuantities.h"
#include "LoKI-B/WorkingConditions.h"
#include "LoKI-B/json.h"

#include <vector>

namespace loki
{

class EedfGasMixture : public GasMixture<EedfGas>
{
  public:
    using Collision = EedfCollision;

    /** Initializes the gas mixture by loading the desired collisions from LXCat files.
     *  These files are read from the electron kinetics setup structure. It also
     *  requires a pointer to the energy grid in order to properly initialize the
     *  cross sections of the collisions.
     */
    EedfGasMixture(Grid *grid, const ElectronKineticsSetup &setup, const WorkingConditions *workingConditions);
    EedfGasMixture(Grid *grid, const json_type &cnf, const WorkingConditions *workingConditions);

    const Vector& elasticCrossSection() const { return m_elasticCrossSection; }
    const Vector& totalCrossSection() const { return m_totalCrossSection; }
    /** \bug The semantics of this member are not clear. Should this consider the 'extra'
     *  collisions as well? At present, they are (in loadCollisions the true flag is set
     *  when a collision of a type is created). At the same time, hasCollisions is used
     *  in ElectronKinetics.cpp to decide whether to add particular terms to the Boltzmann
     *  matrix. If, subsequently, the collsions(that_type) vector is used to populate those
     *  matrix contributions, zero-valued terms may be added (if only *extra* collisions
     *  of that type were loaded).
     */
    bool hasCollisions(CollisionType type) const { return m_hasCollisions[static_cast<uint8_t>(type)]; }
    /// \todo comment evaluateTotalAndElasticCS
    void evaluateTotalAndElasticCS();
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
    void evaluateRateCoefficients(const Vector &eedf);

    std::vector<const EedfGas *> CARGases;
    std::vector<RateCoefficient> rateCoefficients;
    std::vector<RateCoefficient> rateCoefficientsExtra;
  private:
    /** \todo Update docs, this now does ALL the work...
     *  This member contains the bits of the collision creation code
     *  that are shared by the legacy and JSON set up code. It adds
     *  the states that are mentioned on the LHS and RHS of the equation,
     *  then creates the collision object. If an eqiuvalent object
     *  already exists (same particles on both sides, same type), the
     *  object is discarded and a nullptr is returned. Otherwise,
     *  addCollision is called on the target gas of this collision, the
     *  collision is added to our own list of all collisions, the
     *  hasCollisions entry for the particular type of collision is
     *  set to true, and the collision pointer is returned. In this
     *  case, the caller will still need to configure a cross section
     *  object for this collision. That task is not part of this
     *  function since it depends on the input style legacy/JSON.
     */
    void createCollision(const json_type &pcnf, Grid *energyGrid, bool isExtra);

    /** Loads the collisions from the \a file that is provided as first argument.
     *  It also needs a pointer to the \a energyGrid and a boolean \a isExtra to
     *  indicate whether the collisions are extra, for correct initialization and
     *  storage of the collisions.
     *  \todo Explain isExtra.
     */
    void loadCollisions(const std::string &file, Grid *energyGrid, bool isExtra);

    /** Loads the collisions from files, supplied through a vector of strings that hold
     *  the filenames. Furthermore, it needs a pointer to the energy grid and a boolean to
     *  indicate whether the collisions are extra, for correct initialization and storage of
     *  the collisions. When the file extension is ".json", a JSON object is created from
     *  the file and the handling of this file is delegated to member loadCollisionsJSON,
     *  for other file types, the legacy LXCat file format is assumed and member
     *  loadCollisions is called for the file.
     */
    void loadCollisions(const std::vector<std::string> &files, Grid *energyGrid, bool isExtra = false);

    /** EedfGas introduces one extra property that has to be set from a file: OPBParameter.
     *  This override sets this parameter and then calls Gas::loadGasProperties to set the
     *  rest of its properties.
     */
    void loadGasProperties(const GasPropertiesSetup &setup) override;
    /** EedfGas introduces one extra property that has to be set from a file: OPBParameter.
     *  This override sets this parameter and then calls Gas::loadGasProperties to set the
     *  rest of its properties.
     */
    void loadGasProperties(const json_type &cnf) override;

    /// \todo comment addCARGas
    void addCARGas(const std::string& gasName);

    const Grid *grid;

    std::vector<const Collision *> m_collisions;
    Vector m_elasticCrossSection;
    Vector m_totalCrossSection;
    /** \todo Should hasCollisions be recalculated when uMax() is changed?
     *  (smartGrid)? In that case process types may coma and go if they
     *  have threshold above uMax().
     */
    bool m_hasCollisions[static_cast<uint8_t>(CollisionType::size)]{false};
};

} // namespace loki

#endif // LOKI_CPP_EEDFGASMIXTURE_H
