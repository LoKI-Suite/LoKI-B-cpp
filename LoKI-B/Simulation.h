//
// Created by daan on 13-5-19.
//

#ifndef LOKI_CPP_SIMULATION_H
#define LOKI_CPP_SIMULATION_H

#include "LoKI-B/ElectronKinetics.h"
#include "LoKI-B/JobSystem.h"
#include "LoKI-B/Output.h"
#include "LoKI-B/Setup.h"
#include "LoKI-B/WorkingConditions.h"
#include "LoKI-B/json.h"

#include <memory>

namespace loki
{

class Simulation
{
public:
    WorkingConditions m_workingConditions;
    JobManager m_jobManager;
  public:
    ResultEvent m_obtainedResults;
  public:
    explicit Simulation(const Setup &setup);
    Simulation(const json_type &cnf);
    // Copying this object is not allowed.
    Simulation(const Simulation &other) = delete;
    ~Simulation();
    /** Configure an output object. The Simulation object takes ownerwhip
     *  of the pointer.
     */
    void configureOutput(Output* output);

    /// \todo document run
    void run();

  private:
    std::unique_ptr<ElectronKinetics> m_electronKinetics;
    std::unique_ptr<Output> m_output;

    void initializeJobs(const WorkingConditionsSetup &setup);
    void initializeJobs(const json_type &cnf);
};
} // namespace loki

#endif // LOKI_CPP_SIMULATION_H
