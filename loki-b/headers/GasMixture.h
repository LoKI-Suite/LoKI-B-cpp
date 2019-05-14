//
// Created by daan on 14-5-19.
//

#ifndef LOKI_CPP_GASMIXTURE_H
#define LOKI_CPP_GASMIXTURE_H

// TODO: Fill Gas, State and Collision structures
// TODO: Implement LXCat file parsing
// TODO: Implement Property file / function parsing

#include <vector>

#include <EedfGas.h>
#include <EedfState.h>
#include <Collision.h>

namespace loki {
    class GasMixture {
        std::vector<EedfGas> gasses;
        std::vector<EedfState> states;
        std::vector<Collision> collisions, extraCollisions;

        uint32_t addGas(const std::string &name);
//        uint32_t addState(const )
    };
}


#endif //LOKI_CPP_GASMIXTURE_H
