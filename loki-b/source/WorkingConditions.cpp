//
// Created by daan on 2-5-19.
//

#include "WorkingConditions.h"
#include "Constant.h"
#include "Parse.h"
#include "Log.h"
#include "PropertyFunctions.h"

namespace loki {
    using namespace Enumeration;

    // TODO: Comment on WorkingConditions().

    WorkingConditions::WorkingConditions(const WorkingConditionsSetup &setup, const Enumeration::EedfType &eedfType) :
        gasPressure(setup.gasPressure), gasTemperature(setup.gasTemperature),
        electronDensity(setup.electronDensity), chamberLength(setup.chamberLength),
        chamberRadius(setup.chamberRadius), excitationFrequency(setup.excitationFrequency) {

        if ((eedfType == EedfType::boltzmann) &&
                !Parse::getFirstValue(setup.reducedField, reducedField)) {

            Log<ParseFieldError>::Error("reducedField");
        }
        if ((eedfType == EedfType::prescribed) &&
                !Parse::getFirstValue(setup.electronTemperature, electronTemperature)) {

            Log<ParseFieldError>::Error("electronTemperature");
        }

        gasDensity = gasPressure / (Constant::boltzmann * gasTemperature);
        reducedFieldSI = reducedField * 1.e-21;
        reducedExcFreqSI = excitationFrequency * 2 * Constant::pi / gasDensity;

        linkToArgumentMap();
    }

    void WorkingConditions::linkToArgumentMap() {
        argumentMap.emplace("gasTemperature", &gasTemperature);
    }
}