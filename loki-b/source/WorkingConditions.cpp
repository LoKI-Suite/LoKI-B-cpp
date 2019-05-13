//
// Created by daan on 2-5-19.
//

#include <WorkingConditions.h>
#include <Constant.h>

namespace loki {
    // TODO: parse the electronTemperature and reducedField
    //  - We can either parse them in the constructor
    //    (problem is that the job system also needs this information)
    //  - We can initialize the working conditions after creating the jobs.
    //  - We can let the variables be initialized by processing the first job.

    WorkingConditions::WorkingConditions(const WorkingConditionsSetup &setup) :
        gasPressure(setup.gasPressure), gasTemperature(setup.gasTemperature),
        electronDensity(setup.electronDensity), chamberLength(setup.chamberLength),
        chamberRadius(setup.chamberRadius), excitationFrequency(setup.excitationFrequency) {

        gasDensity = gasPressure / (Constant::boltzmann * gasTemperature);
        reducedFieldSI = reducedField * 1.e-21;
        reducedExcFreqSI = excitationFrequency * 2 * Constant::pi / gasDensity;
    }
}