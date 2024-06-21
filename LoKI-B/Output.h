/** \file
 *
 *  Declarations of classes for output generation.
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
 *  \author Jan van Dijk, Daan Boer and Jop Hendrikx
 *  \date   4 July 2019
 */

#ifndef LOKI_CPP_OUTPUT_H
#define LOKI_CPP_OUTPUT_H

#include <string>
#include <functional>

#include "LoKI-B/Grid.h"
#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/MacroscopicQuantities.h"
#include "LoKI-B/Power.h"
#include "LoKI-B/WorkingConditions.h"
#include "LoKI-B/json.h"
#include "LoKI-B/Exports.h"

#ifdef LOKIB_USE_HIGHFIVE

#include <highfive/H5Easy.hpp>

#endif // LOKIB_USE_HDF5

namespace loki
{

class EedfCollisionDataMixture;

class lokib_export Output
{
public:
    virtual ~Output();
    void saveCycle(const Grid &energyGrid, const Vector &eedf, const WorkingConditions &wc, const Power &power,
                   const EedfCollisionDataMixture& collData, const SwarmParameters &swarmParameters,
                   const Vector *firstAnisotropy);
protected:
    Output(const json_type &cnf, const WorkingConditions *workingConditions);
    virtual void setDestination(const std::string& subFolder)=0;
    virtual void writeEedf(const Vector &eedf, const Vector *firstAnisotropy, const Vector &energies) const=0;
    virtual void writeSwarm(const SwarmParameters &swarmParameters) const=0;
    virtual void writePower(const Power &power, const EedfCollisionDataMixture &collData) const=0;
    virtual void writeRateCoefficients(const std::vector<RateCoefficient> &rateCoefficients,
                               const std::vector<RateCoefficient> &extraRateCoefficients) const=0;
    virtual void writeLookupTable(const Power &power,
                                  const std::vector<RateCoefficient> &rateCoefficients,
                                  const std::vector<RateCoefficient> &extraRateCoefficients,
                                  const SwarmParameters &swarmParameters) const=0;
    const WorkingConditions *m_workingConditions;
    bool isBoltzmann() const { return m_isBoltzmann; }
    bool isSimulationHF() const { return m_isSimulationHF; }
private:
    const bool m_isBoltzmann;
    const bool m_isSimulationHF;
    bool saveEedf, savePower, saveSwarm, saveRates, saveTable;
};

class lokib_export FileOutput : public Output
{
public:
    using PathExistsHandler = std::function<void(std::string&)>;
    FileOutput(const json_type &cnf, const WorkingConditions *workingConditions,
        const PathExistsHandler& handler);
protected:
    virtual void setDestination(const std::string& subFolder);
    virtual void writeEedf(const Vector &eedf, const Vector *firstAnisotropy, const Vector &energies) const;
    virtual void writeSwarm(const SwarmParameters &swarmParameters) const;
    virtual void writePower(const Power &power, const EedfCollisionDataMixture &collData) const;
    virtual void writeRateCoefficients(const std::vector<RateCoefficient> &rateCoefficients,
                               const std::vector<RateCoefficient> &extraRateCoefficients) const;
    virtual void writeLookupTable(const Power &power,
                                  const std::vector<RateCoefficient> &rateCoefficients,
                                  const std::vector<RateCoefficient> &extraRateCoefficients,
                                  const SwarmParameters &swarmParameters) const;
private:
    // called by writeLookupTable
    void writeLookupTablePower(const Power &power) const;
    void writeLookupTableRC(const std::vector<RateCoefficient> &rateCoefficients,
                            const std::vector<RateCoefficient> &extraRateCoefficients) const;
    void writeLookupTableSwarmParams(const SwarmParameters &swarmParameters,
                                     const Power &power) const;
    void createPath(const PathExistsHandler& handler);
    void writeTerm(std::ostream& os, const std::string& name, const std::string& unit, double value, bool plus=false) const;
    std::string m_folder;
    std::string m_subFolder;
    mutable bool m_initTable;
};

class lokib_export JsonOutput : public Output
{
public:
    JsonOutput(json_type& root, const json_type &cnf, const WorkingConditions *workingConditions);
protected:
    virtual void setDestination(const std::string& subFolder);
    virtual void writeEedf(const Vector &eedf, const Vector *firstAnisotropy, const Vector &energies) const;
    virtual void writeSwarm(const SwarmParameters &swarmParameters) const;
    virtual void writePower(const Power &power, const EedfCollisionDataMixture& collData) const;
    virtual void writeRateCoefficients(const std::vector<RateCoefficient> &rateCoefficients,
                               const std::vector<RateCoefficient> &extraRateCoefficients) const;
    virtual void writeLookupTable(const Power &power,
                                  const std::vector<RateCoefficient> &rateCoefficients,
                                  const std::vector<RateCoefficient> &extraRateCoefficients,
                                  const SwarmParameters &swarmParameters) const;
private:
    /** This produces a member of the form: name: { "value": value, "unit": unit }.
     *  As an example, makeQuantity("Te", 2.0, "eV") produces and returns an object that contains
     *  \verbatim
          "Te": { "value": = 2.0, "unit": "eV" }. \endverbatim
     *
     */
    static json_type makeQuantity(const std::string& name, double value, const std::string unit);
    json_type& m_root;
    json_type* m_active;
};

#ifdef LOKIB_USE_HIGHFIVE

class lokib_export HDF5Output : public Output
{
public:
    HDF5Output(const std::filesystem::path& fileName, const json_type &cnf, const WorkingConditions *workingConditions);
protected:
    virtual void setDestination(const std::string& subFolder);
    virtual void writeEedf(const Vector &eedf, const Vector *firstAnisotropy, const Vector &energies) const;
    virtual void writeSwarm(const SwarmParameters &swarmParameters) const;
    virtual void writePower(const Power &power, const EedfCollisionDataMixture& collData) const;
    virtual void writeRateCoefficients(const std::vector<RateCoefficient> &rateCoefficients,
                               const std::vector<RateCoefficient> &extraRateCoefficients) const;
    virtual void writeLookupTable(const Power &power,
                                  const std::vector<RateCoefficient> &rateCoefficients,
                                  const std::vector<RateCoefficient> &extraRateCoefficients,
                                  const SwarmParameters &swarmParameters) const;
private:
    /// \todo needed?
    static json_type makeQuantity(const std::string& name, double value, const std::string unit);
    mutable H5Easy::File m_file;
    std::string m_active_path;
};

#endif // LOKIB_USE_HIGHFIVE

} // namespace loki

#endif // LOKI_CPP_OUTPUT_H
