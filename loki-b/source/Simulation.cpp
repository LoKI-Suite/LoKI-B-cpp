//
// Created by daan on 13-5-19.
//

#include "Simulation.h"

// TODO: It might deem necessary to pass the full ElectronKineticsSetup structure
//  to the WorkingConditions constructor (such that we can also check whether it
//  is enabled).

loki::Simulation::Simulation(const loki::Setup &setup)
    : workingConditions(setup.workingConditions, setup.electronKinetics.eedfType),
      enableKinetics(setup.electronKinetics.isOn) {

    if (enableKinetics) {
        electronKinetics = std::make_unique<ElectronKinetics>(setup.electronKinetics, &workingConditions);
    }
}
