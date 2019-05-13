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

#include <Event.h>
#include <Setup.h>

namespace loki {
    class WorkingConditions {
    public:
        double gasPressure,
               gasTemperature,
               gasDensity,
               electronDensity,
               electronTemperature,
               chamberLength,
               chamberRadius,
               reducedField,
               reducedFieldSI,
               excitationFrequency,
               reducedExcFreqSI;

        explicit WorkingConditions(const WorkingConditionsSetup &setup);
        ~WorkingConditions() = default;

        // Copying this object is not allowed.
        WorkingConditions(const WorkingConditions &other) = delete;

        /*
         * We only need one 'updatedGasTemperature' event, since we can control the order
         * of the listeners.
         */

        // Events
        event updatedGasPressure,
              updatedGasTemperature,
              updatedGasDensity,
              updatedElectronDensity,
              updatedElectronTemperature,
              updatedChamberLength,
              updatedReducedField,
              updatedExcitationFrequency;

        // TODO: add functions (or a single templated function of some sort) to update
        //  any variables that can be a range. E.g. at this moment this is 'electronTemperature'
        //  and 'reducedField'.
    };
}


#endif //LOKI_CPP_WORKINGCONDITIONS_H
