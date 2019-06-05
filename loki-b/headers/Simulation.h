//
// Created by daan on 13-5-19.
//

#ifndef LOKI_CPP_SIMULATION_H
#define LOKI_CPP_SIMULATION_H

#include "Setup.h"
#include "WorkingConditions.h"
#include "ElectronKinetics.h"

#include <memory>

namespace loki {
    class Simulation {
        WorkingConditions workingConditions;
        std::unique_ptr<ElectronKinetics> electronKinetics;
        const bool enableKinetics;

    public:
        explicit Simulation(const Setup &setup);
        ~Simulation() = default;

        // Copying this object is not allowed.
        Simulation(const Simulation &other) = delete;

    };
}


#endif //LOKI_CPP_SIMULATION_H
