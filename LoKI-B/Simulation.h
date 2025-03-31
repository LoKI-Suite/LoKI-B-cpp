/** \file
 *
 *  Declaration of LoKI-B's Simulation class.
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
 *  \author Daan Boer and Jan van Dijk (C++ version)
 *  \date   13 May 2019 (first C++ version)
 */


#ifndef LOKI_CPP_SIMULATION_H
#define LOKI_CPP_SIMULATION_H

#include "LoKI-B/ElectronKinetics.h"
#include "LoKI-B/JobSystem.h"
#include "LoKI-B/WorkingConditions.h"
#include "LoKI-B/json.h"
#include "LoKI-B/Exports.h"

#include <memory>

namespace loki
{

/** The simulation class manages all that is needed for a series of runs.
 *
 *  The configuration is described by two data members. The physical properties
 *  of the mixture and the algorithm for evaluating the EEDF and derived
 *  parameters are managed by member m_electronKinetics. The parameters that
 *  describe a particular model run are managed by member m_workingConditions.
 *  That class manages the actual values of parameters such as the actual E/N
 *  value and provides a member for updating that value and executing the tasks
 *  that must be carried out in preparation of a new model run as a result
 *  (WoekingConditions::updateReducedField, in the case of E/N).
 *
 *  Running simulations for one or more parameter ranges is controlled by the
 *  m_jobmanager. It allows the addition of parameters and for each parameter
 *  it records the function that must be invoked when that paramer changes. As
 *  an example, a range of E/N values can be specified together with a pointer
 *  to the updateReducedField member of the WorkingCOnditions object. When Run
 *  is invoked, for each job (combination of parameters) the working conditions
 *  are updated (and the necessary preparatoru tasks executed) and Solve is
 *  called on the electron kinetics object.
 *
 *  Finally, the class provides member m_obtainedResults. In the constructor
 *  of this class, the electron kinetics member is told that the actions that
 *  are registered with this member must be executed when a new run has
 *  completed. After creating a Simulation objec (but before calling Run),
 *  such tasks can be registered (outut tasks, typically). A code may register
 *  multiple handlers with the m_obtainedResults event handler, for example
 *  one for writing/updating on-disk output files, and one for displaying
 *  diagnostic information on the screen. A gui application may add a handler
 *  for updating a plot that is shown on the screen.
 */
class lokib_export Simulation
{
public:
    Simulation(const std::filesystem::path &basePath, const json_type &cnf);
    Simulation(const Simulation &other) = delete;
    Simulation(Simulation &&other) = delete;
    Simulation& operator=(const Simulation &other) = delete;
    Simulation& operator=(Simulation &&other) = delete;
    ~Simulation() = default;
    /// \todo document run
    void run();
    const WorkingConditions& workingConditions() const { return m_workingConditions; }
    ElectronKinetics::ResultEvent& obtainedResults() { return m_obtainedResults; }
private:
    void initializeJobs(const json_type &cnf, bool useReducedFieldParameter);
    WorkingConditions m_workingConditions;
    JobManager m_jobManager;
    std::unique_ptr<ElectronKinetics> m_electronKinetics;
    ElectronKinetics::ResultEvent m_obtainedResults;
};
} // namespace loki

#endif // LOKI_CPP_SIMULATION_H
