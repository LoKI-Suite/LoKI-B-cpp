//
// Created by daan on 13-5-19.
//

#ifndef LOKI_CPP_SIMULATION_H
#define LOKI_CPP_SIMULATION_H

#include "LoKI-B/ElectronKinetics.h"
#include "LoKI-B/JobSystem.h"
#include "LoKI-B/Setup.h"
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
    explicit Simulation(const Setup &setup);
    explicit Simulation(const json_type &cnf);
    // Copying this object is not allowed.
    Simulation(const Simulation &other) = delete;
    ~Simulation();

    /// \todo document run
    void run();

private:

    void initializeJobs(const WorkingConditionsSetup &setup);
    void initializeJobs(const json_type &cnf);

    std::unique_ptr<ElectronKinetics> m_electronKinetics;

public:
    /** \todo Make these private, use read-only accessor where necessary
     *  That avoids accidental modifications of these members that could
     *  break the proper operation of this class.
     */
    WorkingConditions m_workingConditions;
    JobManager m_jobManager;
    ResultEvent m_obtainedResults;
};
} // namespace loki

#endif // LOKI_CPP_SIMULATION_H
