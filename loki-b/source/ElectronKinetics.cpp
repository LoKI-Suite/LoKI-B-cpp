//
// Created by daan on 13-5-19.
//

#include "ElectronKinetics.h"

namespace loki {
    ElectronKinetics::ElectronKinetics(const loki::ElectronKineticsSetup &setup)
        : grid(setup.numerics.energyGrid) {

        mixture.initialize(setup, &grid);

        this->eedfType = setup.eedfType;
        this->shapeParameter = setup.shapeParameter;
        this->ionizationOperatorType = setup.ionizationOperatorType;
        this->growthModelType = setup.growthModelType;
        this->includeEECollisions = setup.includeEECollisions;
    }
}