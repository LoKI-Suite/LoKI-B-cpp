//
// Created by daan on 13-5-19.
//

#ifndef LOKI_CPP_SIMULATION_H
#define LOKI_CPP_SIMULATION_H

#include "Setup.h"
#include "WorkingConditions.h"
#include "ElectronKinetics.h"
#include "JobSystem.h"
#include "Output.h"

#include <memory>

namespace loki {

    class Simulation {
        WorkingConditions workingConditions;
        std::unique_ptr<ElectronKinetics> electronKinetics;
        Output *output;
        const bool enableKinetics, enableOutput;

        bool multipleSimulations{false};

        JobManager jobManager;

    public:
        ResultEvent obtainedResults;
        Event<std::string> outputPathExists;
    public:
        explicit Simulation(const Setup &setup);

        ~Simulation();

        // Copying this object is not allowed.
        Simulation(const Simulation &other) = delete;

        // TODO: comment run

        void run();

    private:
        void initializeJobs(const WorkingConditionsSetup &setup);

        bool initializeJob(const std::string &name, const std::string &valueString,
                           void (WorkingConditions::*callback)(double));
    };
}


#endif //LOKI_CPP_SIMULATION_H
