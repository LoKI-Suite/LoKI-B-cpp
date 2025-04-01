/** \file
 *
 *  Declaration of a class that manages the working conditions of a simulation run.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2025 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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
 *  \author Daan Boer and Jan van Dijk
 *  \date   2 May 2019
 */

#ifndef LOKI_CPP_WORKINGCONDITIONS_H
#define LOKI_CPP_WORKINGCONDITIONS_H

#include "LoKI-B/Event.h"
#include "LoKI-B/json.h"
#include "LoKI-B/Constant.h"
#include <string>
#include <map>

namespace loki
{

/** The WorkingConditions class stores all global variables
 *  that make up the working environment of the simulation.
 *  It is also responsible for propagating any change in the
 *  working conditions to the relevant classes.
 */
class WorkingConditions
{
  public:

    // constructors and destructor

    WorkingConditions(const json_type &cnf);
    ~WorkingConditions() = default;
    // Copying this object is not allowed.
    WorkingConditions(const WorkingConditions &other) = delete;

    /// \todo This is input only for the prescribed EEDF case
    double electronTemperature() const { return m_electronTemperature; }
    double gasTemperature() const { return m_gasTemperature; }
    double gasDensity() const { return m_gasDensity; }
    double electronDensity() const { return m_electronDensity; }
    double reducedField() const { return m_reducedField; }
    double reducedFieldSI() const { return m_reducedField * SI::Townsend; }
    double excitationFrequency() const { return m_excitationFrequency; }
    double reducedExcFreqSI() const { return m_excitationFrequency * 2 * Constant::pi / m_gasDensity; }

    const double* findParameter(const std::string& name) const;
    const double& getParameter(const std::string& name) const;

    void updateReducedField(double value);
    void updateElectronTemperature(double value);
    /** Return a string-identifier that identifies the present case (set of
     *  working conditions, that is). This can be used in output generation
     *  tasks. As an example, file output can use the identifier as a directory
     *  name in which the results for this run are stored, while JSON output can
     *  use the identifier as the nama of the section for that output.
     *
     *  Note that at present the m_currentJobFolder is not updated automatically
     *  by this class when, for example, a new value of the reduced field is
     *  configured. It must rather be set by client code after new parameters
     *  have been configured.
     *
     *  \sa setCurrentJobFolder
     */
    std::string getCurrentJobFolder() const { return m_currentJobFolder; }
    /** Set the string that identifies the present case. See member
     *  setCurrentJobFolder for further details.
     *
     *  \sa getCurrentJobFolder
     */
    void setCurrentJobFolder(const std::string& folder) { m_currentJobFolder = folder; }

    /** \todo Technically we only need one 'updatedGasTemperature' event,
     *        since we can control the order of the listeners.
     */
    /** \todo Document the individual events carefully. What are the intent
     *        and expectations?
     */
    /** \todo Check which events are actually used. Note that this class has
     *        only two 'update' members, for the reduced field and for the
     *        electron temperature. Update: only those two events are used
     *        at present.
     */
    Event<> updatedReducedField;
    Event<> updatedElectronTemperature;
    //Event<> updatedGasPressure;
    //Event<> updatedGasTemperature1;
    //Event<> updatedGasTemperature2;
    //Event<> updatedGasDensity;
    //Event<> updatedElectronDensity;
    //Event<> updatedChamberLength;
    //Event<> updatedExcitationFrequency;
private:
    const double* addParameter(const std::string& name, const double* vptr);

    const double m_gasPressure;
    const double m_gasTemperature;
    const double m_gasDensity;
    double m_electronDensity;
    /** \todo It appears that electronTemperature is not used anywhere
     *  at present, maybe because only EEDF type 'boltzmann' has been
     *  implemented? It is not used, and set only by
     *  evaluateSwarmParameters() at present. Check this.
     */
    double m_electronTemperature;
    /// \todo chamberLength is not used anywhere at present
    //double chamberLength;
    /// \todo chamberRadius is not used anywhere at present
    //double chamberRadius;
    double m_reducedField;
    double m_excitationFrequency;
    using ArgumentMap = std::map<std::string, const double *>;
    ArgumentMap m_argumentMap;
    std::string m_currentJobFolder;
};

} // namespace loki

#endif // LOKI_CPP_WORKINGCONDITIONS_H
