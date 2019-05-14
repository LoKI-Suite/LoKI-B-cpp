//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_GAS_H
#define LOKI_CPP_GAS_H

#include <string>

namespace loki {
    class Gas {
        std::string name;
        double mass,
               harmonicFrequency,
               anharmonicFrequency,
               rotationalConstant,
               electricDipoleMoment,
               electricQuadrupoleMoment,
               polarizability,
               fraction;

    public:
        bool operator==(const Gas &other);
    };
}


#endif //LOKI_CPP_GAS_H
