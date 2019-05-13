//
// Created by daan on 13-5-19.
//

#ifndef LOKI_CPP_ELECTRONKINETICS_H
#define LOKI_CPP_ELECTRONKINETICS_H

#include <Setup.h>

namespace loki {
    class ElectronKinetics {

    public:
        explicit ElectronKinetics(const ElectronKineticsSetup &setup);
        ~ElectronKinetics() = default;

        // Copying this object is not allowed.
        ElectronKinetics(const ElectronKinetics &other) = delete;
    };
}


#endif //LOKI_CPP_ELECTRONKINETICS_H
