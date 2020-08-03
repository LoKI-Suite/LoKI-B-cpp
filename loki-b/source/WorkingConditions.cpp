//
// Created by daan on 2-5-19.
//

#include "LoKI-B/WorkingConditions.h"
#include "LoKI-B/Constant.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/PropertyFunctions.h"
#include <limits>

namespace loki {
    using namespace Enumeration;

    // TODO: Comment on WorkingConditions().

    WorkingConditions::WorkingConditions(const WorkingConditionsSetup &setup) :
        gasPressure(setup.gasPressure),
        gasTemperature(setup.gasTemperature),
        electronDensity(setup.electronDensity),
        chamberLength(setup.chamberLength),
        chamberRadius(setup.chamberRadius),
        excitationFrequency(setup.excitationFrequency)
    {
        /* set the reducedField and electronTemperature to dummy
         * values. These are set by the JobManager when prepareFirsJobs
         * or nextJov() is called.
         */
        reducedField = std::numeric_limits<double>::quiet_NaN();
        /** \todo It appears that electronTemperature is not used anywhere
         *  at present, maybe because only EEDF type 'boltzmann' has been
         *  implemented? It is not used, and set only by
         *  evaluateSwarmParameters() at present. Check this.
         */
        electronTemperature = std::numeric_limits<double>::quiet_NaN();

        gasDensity = gasPressure / (Constant::boltzmann * gasTemperature);
        reducedFieldSI = reducedField * 1.e-21;
        reducedExcFreqSI = excitationFrequency * 2 * Constant::pi / gasDensity;

        linkToArgumentMap();
    }

    WorkingConditions::WorkingConditions(const json_type &cnf) :
        gasPressure(cnf.at("gasPressure").get<double>()),
        gasTemperature(cnf.at("gasTemperature").get<double>()),
        electronDensity(cnf.at("electronDensity").get<double>()),
        chamberLength(cnf.at("chamberLength").get<double>()),
        chamberRadius(cnf.at("chamberRadius").get<double>()),
        excitationFrequency(cnf.at("excitationFrequency").get<double>()) {

        /* set the reducedField and electronTemperature to dummy
         * values. These are set by the JobManager when prepareFirsJobs
         * or nextJov() is called.
         */
        reducedField = std::numeric_limits<double>::quiet_NaN();
        /** \todo It appears that electronTemperature is not used anywhere
         *  at present, maybe because only EEDF type 'boltzmann' has been
         *  implemented? It is not used, and set only by
         *  evaluateSwarmParameters() at present. Check this.
         */
        electronTemperature = std::numeric_limits<double>::quiet_NaN();

        gasDensity = gasPressure / (Constant::boltzmann * gasTemperature);
        reducedFieldSI = reducedField * 1.e-21;
        reducedExcFreqSI = excitationFrequency * 2 * Constant::pi / gasDensity;

        linkToArgumentMap();
    }

    void WorkingConditions::updateReducedField(double value) {
        reducedField = value;
        reducedFieldSI = reducedField * 1.e-21;

        Log<Message>::Notify(value);

        updatedReducedField.emit();
    }

    void WorkingConditions::updateElectronTemperature(double value) {
        electronTemperature = value;

        updatedElectronTemperature.emit();
    }

    void WorkingConditions::linkToArgumentMap() {
        argumentMap.emplace("gasTemperature", &gasTemperature);
    }
}
