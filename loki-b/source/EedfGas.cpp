//
// Created by daan on 2-5-19.
//

#include <EedfGas.h>

namespace loki {
    EedfGas::EedfGas(const std::string &name) : Gas(name) {}

    void EedfGas::addCollision(EedfCollision *collision, bool isExtra) {

        (isExtra ? this->extraCollisions : this->collisions).emplace_back(collision);
    }

    EedfGas::~EedfGas() {
        for (auto *collision : collisions)
            delete collision;

        for (auto *collision : extraCollisions)
            delete collision;
    }
}