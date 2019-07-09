//
// Created by daan on 13-5-19.
//

#include "Simulation.h"
#include <chrono>

namespace loki {
    Simulation::Simulation(const loki::Setup &setup)
            : workingConditions(setup.workingConditions, setup.electronKinetics.eedfType),
              enableKinetics(setup.electronKinetics.isOn), enableOutput(setup.output.isOn),
              jobManager(&workingConditions) {

        if (enableKinetics) {
            initializeJobs(setup.workingConditions);

            electronKinetics = std::make_unique<ElectronKinetics>(setup.electronKinetics, &workingConditions);
            electronKinetics->obtainedNewEedf.addListener(&ResultEvent::emit, &obtainedResults);

            if (enableOutput) {
                output = new Output(setup.output, &workingConditions, &jobManager);

                electronKinetics->obtainedNewEedf.addListener(&Output::saveCycle, output);
                output->simPathExists.addListener(&Event<std::string>::emit, &outputPathExists);
            }
        }
    }

    void Simulation::run() {
        if (enableKinetics) {
            if (enableOutput) {
                output->createPath();
            }
            if (multipleSimulations) {
                do {
                    electronKinetics->solve();
                } while (jobManager.nextJob());
            } else {
                electronKinetics->solve();
            }

        }
    }

    Simulation::~Simulation() {
        if (enableOutput) delete output;
    }

    void Simulation::initializeJobs(const WorkingConditionsSetup &setup) {

        // Repeat this statement for any other fields that can be declared as a range.
        if (!Parse::isNumerical(setup.reducedField)) {
            if (!initializeJob("Reduced Field", setup.reducedField, &WorkingConditions::updateReducedField)) {
                Log<Message>::Error("Reduced field entry in input file is ill formatted.");
            }
        }
    }

    bool Simulation::initializeJob(const std::string &name, const std::string &valueString,
                                   void (WorkingConditions::*callback)(double)) {
        bool success = false;

        Range range = Parse::getRange(valueString, success);

        if (!success) return false;

        jobManager.addJob({name, &WorkingConditions::updateReducedField, range});
        multipleSimulations = true;

        return true;
    }
}