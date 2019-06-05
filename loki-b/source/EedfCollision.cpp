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
            Log<MultipleReactantInEedfCol>::Warning(*this);

        if (isReverse && products.size() > 1)
            Log<MultipleProductsInReverse>::Error(*this);
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
           << (collision.isReverse ? " <->" : " ->");

        if (collision.type != CollisionType::attachment) os << " e +";
        if (collision.type == CollisionType::ionization) os << " e +";

        for (uint32_t i = 0; i < collision.products.size(); ++i) {
            os << ' ' << *collision.products[i];
            if (i < collision.products.size() - 1) os << " +";
        }

        return os;
    }

    void EedfCollision::superElastic(const Vector &energyData, Vector &result) const {
        if (!isReverse) {
            Log<SuperElasticForNonReverse>::Error(*this);
        }

        result.resize(energyData.size());

        Vector superElasticEnergies = energyData.array() + crossSection->threshold;

        crossSection->interpolate(superElasticEnergies, result);

        const double swRatio = target->statisticalWeight / products[0]->statisticalWeight;

        uint32_t minIndex = 0;

        if (energyData[0] == 0)
            ++minIndex;

        for (uint32_t i = minIndex; i < result.size(); ++i) {
            result[i] *= swRatio * (1 + (crossSection->threshold / energyData[i]));
        }

        if (energyData[0] == 0) result[0] = 0;
    }
}