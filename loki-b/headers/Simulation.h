//
// Created by daan on 13-5-19.
//

#ifndef LOKI_CPP_SIMULATION_H
#define LOKI_CPP_SIMULATION_H

#include "Setup.h"
#include "json.h"
#include "WorkingConditions.h"
#include "ElectronKinetics.h"
#include "JobSystem.h"
#include "Output.h"

#include <memory>

namespace loki {

    class Simulation {
        WorkingConditions workingConditions;
        std::unique_ptr<ElectronKinetics> electronKinetics;
        std::unique_ptr<Output> output;
        JobManager jobManager;
    public:
        ResultEvent obtainedResults;
        Event<std::string> outputPathExists;
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
    };
}


#endif //LOKI_CPP_SIMULATION_H
