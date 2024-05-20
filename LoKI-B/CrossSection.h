/** \file
 *
 *  Declaration of a class that represents cross sections.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2024 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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
 *  \date   15 May 2019
 */

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
    CrossSection(const Grid *energyGrid, bool isElasticOrEffective, const json_type &cnf);
    CrossSection(double threshold, const Grid *energyGrid, bool isElasticOrEffective, std::istream &in);
    CrossSection(double threshold, const Grid *energyGrid, bool isElasticOrEffective, Vector rawEnergyData,
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
