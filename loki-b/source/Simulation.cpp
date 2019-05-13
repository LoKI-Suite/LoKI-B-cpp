//
// Created by daan on 13-5-19.
//

#include "Simulation.h"

loki::Simulation::Simulation(const loki::Setup &setup)
    : workingConditions(setup.workingConditions),
      enableKinetics(setup.electronKinetics.isOn) {

    if (enableKinetics) {
        electronKinetics = new ElectronKinetics(setup.electronKinetics);
    }

}

loki::Simulation::~Simulation() {

}
