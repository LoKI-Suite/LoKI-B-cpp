//
// Created by daan on 15-5-19.
//

#ifndef LOKI_CPP_CROSSSECTION_H
#define LOKI_CPP_CROSSSECTION_H

#include "LoKI-B/Grid.h"
#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/json.h"
#include "LoKI-B/LookupTable.h"

#include <iostream>
#include <vector>

namespace loki
{

class CrossSection : public Vector
{
public:
    CrossSection(Grid *energyGrid, bool isElasticOrEffective, const json_type &cnf);
    CrossSection(double threshold, Grid *energyGrid, bool isElasticOrEffective, std::istream &in);
    CrossSection(double threshold, Grid *energyGrid, bool isElasticOrEffective, Vector rawEnergyData,
                 Vector rawCrossSection);

    void interpolate();
    void interpolate(const Vector &energies, Vector &result) const;
    const Grid *getGrid() const { return m_energyGrid; }
    double threshold() const { return m_threshold; }

    const LookupTable& lookupTable() const { return m_lut; }
private:
    using Index = Vector::Index;
    const double m_threshold;
    const Grid *m_energyGrid;
    const bool m_isElasticOrEffective;
    const LookupTable m_lut;
};

} // namespace loki

#endif // LOKI_CPP_CROSSSECTION_H
