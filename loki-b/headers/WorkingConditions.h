//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_WORKINGCONDITIONS_H
#define LOKI_CPP_WORKINGCONDITIONS_H

/*
 * The WorkingConditions class stores all global variables
 * that make up the working environment of the simulation.
 * It is also responsible for propagating any change in the
 * working conditions to the relevant classes.
 */

#include "Event.h"
#include "Setup.h"
#include "Enumeration.h"
#include <map>

namespace loki {
    class WorkingConditions {
    public:
        double gasPressure,
               gasTemperature,
               gasDensity,
               electronDensity,
               electronTemperature{0},
               chamberLength,
               chamberRadius,
               reducedField{0},
               reducedFieldSI,
               excitationFrequency,
               reducedExcFreqSI;

        std::map<std::string, double *> argumentMap;

        explicit WorkingConditions(const WorkingConditionsSetup &setup, const Enumeration::EedfType &eedfType);
        ~WorkingConditions() = default;

        // Copying this object is not allowed.
        WorkingConditions(const WorkingConditions &other) = delete;

        void updateReducedField(double value);

        void updateElectronTemperature(double value);

        /*
         * Technically we only need one 'updatedGasTemperature' event, since we can control
         * the order of the listeners.
         */

        // Events
        event updatedGasPressure,
              updatedGasTemperature1,
              updatedGasTemperature2,
              updatedGasDensity,
              updatedElectronDensity,
              updatedElectronTemperature,
              updatedChamberLength,
              updatedReducedField,
              updatedExcitationFrequency;

    private:
        void linkToArgumentMap();
    };
}


#endif //LOKI_CPP_WORKINGCONDITIONS_H
