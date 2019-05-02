//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_SETUP_H
#define LOKI_CPP_SETUP_H

#include <string>
#include <vector>
#include <Enumerations.h>

namespace loki {
    // TODO: think about whether this is the right thing to do,
    //  would it not be more convenient to store the values
    //  immediately in their designated simulation structures?

    struct WorkingConditionsSetup {
        std::string reducedField;
        std::string electronTemperature;
        double excitationFrequency = 0.;
        double gasPressure = 0.;
        double gasTemperature = 0.;
        double electronDensity = 0.;
        double chamberLength = 0.;
        double chamberRadius = 0.;
    };

    struct electronKineticsSetup {
        bool isEnabled = false;
        EedfType eedf;
        uint8_t shapeParameter = 0;
        IonizationOperatorType ionizationOperator;
        GrowthModelType growthModel;
        bool includeEECollisions = false;
        std::vector<std::string> LXCatFiles;
        std::vector<std::string> extraLXCatFiles;
    };

    struct Setup {

    };


}


#endif //LOKI_CPP_SETUP_H
