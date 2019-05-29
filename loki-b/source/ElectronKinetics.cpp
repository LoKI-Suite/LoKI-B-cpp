//
// Created by daan on 13-5-19.
//

#include "ElectronKinetics.h"
#include <chrono>

namespace loki {
    ElectronKinetics::ElectronKinetics(const ElectronKineticsSetup &setup, const WorkingConditions *workingConditions)
        : workingConditions(workingConditions), grid(setup.numerics.energyGrid) {

        mixture.initialize(setup, &grid, workingConditions);

        this->eedfType = setup.eedfType;
        this->shapeParameter = setup.shapeParameter;
        this->ionizationOperatorType = setup.ionizationOperatorType;
        this->growthModelType = setup.growthModelType;
        this->includeEECollisions = setup.includeEECollisions;
    }
}