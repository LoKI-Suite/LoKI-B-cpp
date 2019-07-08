//
// Created by daan on 13-5-19.
//

#include "Simulation.h"
#include <chrono>

// TODO: It might deem necessary to pass the full ElectronKineticsSetup structure
//  to the WorkingConditions constructor (such that we can also check whether it
//  is enabled).
// TODO [FUTURE]: We might want to extend the range behaviour of the reduced field to other settings
//  in the working conditions. In this case one can link a function to update the corresponding member variable
//  to a range and let a 'manager' class handle the next value to update (such that all combinations of values
//  appear once).

namespace loki {
    Simulation::Simulation(const loki::Setup &setup)
            : workingConditions(setup.workingConditions, setup.electronKinetics.eedfType),
              enableKinetics(setup.electronKinetics.isOn), enableOutput(setup.output.isOn),
              jobManager(&workingConditions) {

        if (enableKinetics) {
//            initializeFieldRange(setup.workingConditions);
            initializeJobs(setup.workingConditions);

            electronKinetics = std::make_unique<ElectronKinetics>(setup.electronKinetics, &workingConditions);

            if (enableOutput) {
                output = new Output(setup.output, electronKinetics->getGrid(), &workingConditions, &jobManager);

                electronKinetics->obtainedNewEedf.addListener(&Output::saveCycle, output);
            }
        }
    }

    void Simulation::run() {
        if (enableKinetics) {
            if (multipleSimulations) {
                do {
                    // DONE: update working conditions with the next value in the range.
                    // TODO: update the output class so it creates subfolders.
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