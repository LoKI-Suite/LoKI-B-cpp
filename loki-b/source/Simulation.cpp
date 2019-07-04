//
// Created by daan on 13-5-19.
//

#include "Simulation.h"
#include <chrono>

// TODO: It might deem necessary to pass the full ElectronKineticsSetup structure
//  to the WorkingConditions constructor (such that we can also check whether it
//  is enabled).

namespace loki {
    Simulation::Simulation(const loki::Setup &setup)
            : workingConditions(setup.workingConditions, setup.electronKinetics.eedfType),
              enableKinetics(setup.electronKinetics.isOn), enableOutput(setup.output.isOn) {

        if (enableKinetics) {
            electronKinetics = std::make_unique<ElectronKinetics>(setup.electronKinetics, &workingConditions);

            if (enableOutput) {
                output = new Output(setup.output, electronKinetics->getGrid(), &workingConditions);

                electronKinetics->obtainedNewEedf.addListener(&Output::saveCycle, output);
            }
        }
    }

    void Simulation::run() {
        if (enableKinetics) {
            electronKinetics->solve();
        }
    }

    Simulation::~Simulation() {
        delete output;
    }
}