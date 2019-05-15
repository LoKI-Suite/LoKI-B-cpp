//
// Created by daan on 13-5-19.
//

#include "ElectronKinetics.h"

namespace loki {
    ElectronKinetics::ElectronKinetics(const loki::ElectronKineticsSetup &setup) : grid(setup.numerics.energyGrid) {}
}