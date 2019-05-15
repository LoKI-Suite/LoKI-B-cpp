//
// Created by daan on 15-5-19.
//

#include <iostream>

#include <utility>
#include <vector>
#include <string>
#include <algorithm>
#include <map>

class Boltzmann;
class Chemistry;

class EedfGas;
class EedfState;

class ChemGas;
class ChemState;

enum StateLevel : uint8_t {
    electronic,
    vibrational,
    rotational,
    ionic,
    nan
};

struct StateEntry {
    StateLevel level;
    std::string gasName, e;
    // v, J, charge
    int16_t levels[3];
};

template <typename TraitType> struct Trait;

template <>
struct Trait<Boltzmann> {
    typedef EedfGas Gas;
    typedef EedfState State;
};

template <>
struct Trait<Chemistry> {
    typedef ChemGas Gas;
    typedef ChemState State;
};

template <typename TraitType>
struct Gas {
    const std::string name;
    std::map<std::string, typename Trait<TraitType>::State> states;

    explicit Gas(std::string name) : name(std::move(name)) {}

    bool operator==(const std::string &otherName) {
        return name == otherName;
    }

    void print() const {
        std::cout << "Gas: " << name << std::endl;
        for (const auto & [key, state] : states) {
            state.print();
        }
    }
};

struct EedfGas : public Gas<Boltzmann> {
    void someExtraFunction() {
        std::cout << "I am an EedfGas." << std::endl;
    }

    explicit EedfGas(const std::string &name)
        : Gas<Boltzmann>(name) {}
};

struct ChemGas : public Gas<Chemistry> {
    void someExtraFunction() {
        std::cout << "I am a ChemGas." << std::endl;
    }

    explicit ChemGas(const std::string &name)
            : Gas<Chemistry>(name) {}
};

template<typename TraitType>
struct State {
    typename Trait<TraitType>::State *parent{nullptr};
    typename Trait<TraitType>::Gas *gas{nullptr};

    StateLevel level;
    std::string e;
    uint16_t v, J;
    int16_t charge;

    std::vector<typename Trait<TraitType>::State *> children;

    explicit State(const StateEntry &entry, StateLevel lvl = nan) : e(entry.e),
        v(entry.levels[0]), J(entry.levels[1]), charge(entry.levels[2]) {

        if (lvl == nan)
            this->level = entry.level;
        else
            this->level = lvl;
    }

    virtual ~State() {
        for (auto *child : children)
            delete child;
    }

    void print() const {
        std::cout << "level: e = " << e;
        if (level >= vibrational) std::cout << " v = " << v;
        if (level >= rotational) std::cout << " J = " << J;
        if (level == ionic) std::cout << " charge = " << charge;
        std::cout << std::endl;

        for (auto state : children) {
            if (state != nullptr) state->print();
        }
    }
};

struct EedfState : public State<Boltzmann> {
    ChemState *chemEquivalent{nullptr};

    void someExtraFunction() {
        std::cout << "I am an EedfState." << std::endl;
    }

    explicit EedfState(const StateEntry &entry, StateLevel level = nan)
                : State<Boltzmann>(entry, level) {}
};

struct ChemState : public State<Chemistry> {
    EedfState *eedfEquivalent{nullptr};

    void someExtraFunction() {
        std::cout << "I am a ChemState." << std::endl;
    }

    explicit ChemState(const StateEntry &entry, StateLevel level = nan)
                : State<Chemistry>(entry, level) {}
};

template <typename TraitType>
struct GasMixture {
    std::vector<typename Trait<TraitType>::Gas> gasses;

    typename Trait<TraitType>::Gas &addGas(const std::string &name) {
        auto it = std::find(gasses.begin(), gasses.end(), name);

        if (it == gasses.end()) {
            return gasses.emplace_back(name);
        }

        return *it;
    }

    typename Trait<TraitType>::State &
    addState(const StateEntry &entry) {
        auto &gas = addGas(entry.gasName);

        auto *state = &addState(gas, entry);

        if (entry.level == electronic) return *state;

        for (uint8_t l = electronic; l < entry.level; ++l) {
            state = addState(*state, entry);
        }

        return *state;
    }

    typename Trait<TraitType>::State &
    addState(typename Trait<TraitType>::Gas &gas, const StateEntry &entry) {
        return (*(gas.states.emplace(std::piecewise_construct, std::forward_as_tuple(entry.e),
                std::forward_as_tuple(entry, electronic)).first)).second;
    }

    typename Trait<TraitType>::State *
    addState(typename Trait<TraitType>::State &state, const StateEntry &entry) {
        auto &l = entry.levels[state.level];

        if (state.children.size() <= l) {
            state.children.resize(l+1);
        }

        if (state.children[l] == nullptr) {
            state.children[l] = new typename Trait<TraitType>::State(entry, (StateLevel)(state.level+1));
        }

        return state.children[l];
    }

    void print() {
        for (int i = 0; i < gasses.size(); ++i) {
            gasses[i].print();
        }
    }
};

int main (int argc, char ** argv) {
    GasMixture<Boltzmann> mixture;

    StateEntry entry{rotational, "N2", "X", 3, 5, 0};
    StateEntry entry2{rotational, "N2", "X", 3, 8, 0};
    StateEntry entry3{ionic, "N2", "A", 9, 3, 1};

    mixture.addState(entry);
    mixture.addState(entry2);
    mixture.addState(entry3);

//    mixture.gasses.at(0).states.at("X").children.at(3)->print();

    mixture.print();

    return 0;
}