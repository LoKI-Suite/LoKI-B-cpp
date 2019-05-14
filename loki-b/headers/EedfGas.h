//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_EEDFGAS_H
#define LOKI_CPP_EEDFGAS_H

#include "Gas.h"
#include "EedfState.h"

#include <vector>

namespace loki {
    class EedfGas : public Gas {
        std::vector<EedfState> states;
    };
}


#endif //LOKI_CPP_EEDFGAS_H
