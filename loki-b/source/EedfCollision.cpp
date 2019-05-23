//
// Created by daan on 21-5-19.
//

#include "EedfCollision.h"
#include "Log.h"

namespace loki {

    EedfCollision::EedfCollision(Enumeration::CollisionType type, std::vector<EedfState *> &reactants,
                                 std::vector<EedfState *> &products, std::vector<uint16_t> &stoiCoeff, bool isReverse)
                                 : Collision(type, reactants[0], products, stoiCoeff, isReverse) {

        if (reactants.size() != 1)
            Log<Message>::Warning("EedfCollision with more than one reactant!");
    }

    EedfCollision::~EedfCollision() {
        delete crossSection;
    }

    bool EedfCollision::operator==(const EedfCollision &other) {
        if (type != other.type) return false;
        if (target != other.target) return false;

        // Comparing pointers since there is only one instantiation of each state.
        for (const auto *state : other.products) {
            if (std::find(products.begin(), products.end(), state) == products.end())
                return false;
        }

        return true;
    }

    EedfState *EedfCollision::getTarget() {
        return target;
    }

    std::ostream &operator<<(std::ostream &os, const EedfCollision &collision) {
        os << "e + " << *collision.target
                     << (collision.isReverse ? " <->":" ->");

        if (collision.type != CollisionType::attachment) os << " e +";
        if (collision.type == CollisionType::ionization) os << " e +";

        for (uint32_t i = 0; i < collision.products.size(); ++i) {
            os << ' ' << *collision.products[i];
            if (i < collision.products.size()-1) os << " +";
        }

        return os;
    }
}