//
// Created by daan on 21-5-19.
//

#ifndef LOKI_CPP_EEDFMIXTURE_H
#define LOKI_CPP_EEDFMIXTURE_H

#include "LoKI-B/GasMixture.h"
#include "LoKI-B/EedfCollisions.h"

#include <string>
#include <vector>

namespace loki
{

class EedfMixture
{
  public:
    /** Initializes the gas mixture by loading the desired collisions from LXCat files.
     *  These files are read from the electron kinetics setup structure \a cnf. It also
     *  requires a pointer to the energy grid in order to properly initialize the
     *  cross sections of the collisions.
     */
    EedfMixture(const std::filesystem::path &basePath, const Grid *grid, const json_type &cnf, const WorkingConditions *workingConditions);

    const GasMixture& composition() const { return m_composition; }

    const std::vector<const Gas *>& CARGases() const { return m_CARGases; }
    const EedfCollisionDataMixture& collision_data() const { return m_collision_data; }
    EedfCollisionDataMixture& collision_data() { return m_collision_data; }
  private:
    /** Loads the collisions from files, supplied through a vector of strings that hold
     *  the filenames. Furthermore, it needs a pointer to the energy grid and a boolean to
     *  indicate whether the collisions are extra, for correct initialization and storage of
     *  the collisions. When the file extension is ".json", a JSON object is created from
     *  the file and the handling of this file is delegated to member loadCollisionsJSON,
     *  otherwise, the legacy LXCat file format is assumed and member loadCollisionsClassic
     *  is called on the file.
     */
    void loadCollisions(const std::filesystem::path &basePath, const std::vector<std::string> &files, const GasProperties& gasProps, const Grid *energyGrid, bool isExtra = false);

    /// \todo comment addCARGas
    void addCARGas(const std::string& gasName);

    GasMixture m_composition;
    EedfCollisionDataMixture m_collision_data;

    std::vector<const Gas *> m_CARGases;
};

} // namespace loki

#endif // LOKI_CPP_EEDFMIXTURE_H
