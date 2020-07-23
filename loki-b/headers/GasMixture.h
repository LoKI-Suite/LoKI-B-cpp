//
// Created by daan on 14-5-19.
//

#ifndef LOKI_CPP_GASMIXTURE_H
#define LOKI_CPP_GASMIXTURE_H

// DONE: Fill the Gas, State and Collision structures.
// DONE: Implement the templated approach to gasses and states
// DONE: Think of how we can parse the LXCat data. Especially for states, do we need an
//  auxiliary structure like StateEntry or can we immediately create the state?
//  If not then alter the GasMixture class (replace StateEntry by State).
// DONE: Implement LXCat file parsing.
// DONE: Implement Property file / function / value parsing.

#include "EedfGas.h"
#include "EedfState.h"
#include "EedfCollision.h"
#include "Collision.h"
#include "Setup.h"
#include "Traits.h"
#include "Log.h"
#include "PropertyFunctions.h"
#include "WorkingConditions.h"
#include "Parse.h"
#include "json.h"

#include <vector>
#include <set>
#include <fstream>
#include <regex>

/* -- GAS_PROPERTY --
 * This define accepts the name of a gas property (e.g. mass) and will try load the corresponding database
 * file. If this fails it will output a warning (but nothing more). If the file is succesfully loaded, it
 * will then try to set the property for all gasses in the mixture. This version will not throw an error
 * in any case, and can therefore be used for gas properties that are not mandatory for every gas (e.g.
 * atomic gasses do not have a harmonicFrequency).
 *
 * It assumes that the passed property has the same name in the Gas and GasPropertiesSetup structures.
 * Furthermore, it assumes that a std::string fileBuffer has been declared beforehand.
 */

#define GAS_PROPERTY_TEMPL(PROPERTY,SEVERITY) \
{ \
    std::string fileBuffer; \
    std::cout << "Configuring gas property '" << #PROPERTY << "', using file '" << setup.PROPERTY << "'." << std::endl; \
    if (Parse::stringBufferFromFile(setup.PROPERTY, fileBuffer)) { \
        for (auto *gas : gasses) { \
            if(!Parse::gasProperty(gas->name, gas->PROPERTY, fileBuffer)) { \
                Log<GasPropertyError>::SEVERITY(#PROPERTY " in gas " + gas->name); \
            } \
        } \
    } \
    else { \
        Log<FileError>::SEVERITY(setup.PROPERTY); \
    } \
}

#define GAS_PROPERTY(PROPERTY) GAS_PROPERTY_TEMPL(PROPERTY,Warning)
#define R_GAS_PROPERTY(PROPERTY) GAS_PROPERTY_TEMPL(PROPERTY,Error)

// JSON versions. The filename is obtained from cnf, the JSON version of the input file
/// \todo Error reporting should be fixed. The file name is not displayed, the field instead.
#define GAS_PROPERTY_JSON_TEMPL(cnf_object,PROPERTY,SEVERITY) \
{ \
    std::string fileBuffer; \
    if (cnf_object.contains(#PROPERTY) && Parse::stringBufferFromFile(cnf_object.at(#PROPERTY), fileBuffer)) { \
        std::cout << "Configuring gas property '" << #PROPERTY << "', using file '" << cnf_object.at(#PROPERTY) << "'." << std::endl; \
        for (auto *gas : gasses) { \
            if(!Parse::gasProperty(gas->name, gas->PROPERTY, fileBuffer)) { \
                Log<GasPropertyError>::SEVERITY(#PROPERTY " in gas " + gas->name); \
            } \
        } \
    } \
    else { \
        Log<FileError>::SEVERITY(#PROPERTY); \
    } \
}

#define GAS_PROPERTY_JSON(cnf_object,PROPERTY) GAS_PROPERTY_JSON_TEMPL(cnf_object,PROPERTY,Warning);
#define R_GAS_PROPERTY_JSON(cnf_object,PROPERTY) GAS_PROPERTY_JSON_TEMPL(cnf_object,PROPERTY,Error);

namespace loki {
    /* -- GasMixture --
     * The GasMixture class acts as a base class to the EedfGasMixture and future
     * ChemGasMixture classes. This class is templated to allow the use of trait
     * classes (classes that define types). Its constructor is protected such
     * that the class cannot be instantiated as is, although other classes can be
     * derived from this class. They can either inherit from GasMixture<Boltzmann>
     * or GasMixture<Chemistry>. The idea is that the first handles EedfGasses,
     * States and Collisions, whereas the latter handles their Chemistry
     * equivalents.
     */

    template<typename TraitType>
    struct GasMixture {
        // Vector of pointers to all the gasses in the mixture.
        std::vector<typename Trait<TraitType>::Gas *> gasses;

        virtual ~GasMixture() {
            for (auto *gas : gasses) {
                delete gas;
            }
        };

        void print() {
            for (int i = 0; i < gasses.size(); ++i) {
                std::cout << "Gas: " << gasses[i]->name << std::endl;
                gasses[i]->print();
            }
        }

    protected:
        GasMixture() = default;

        /* -- findGas --
         * Tries to find a gas in the gas mixture specified by a supplied name. If the
         * gas is present, a pointer to this gas is returned, otherwise a nullptr is
         * returned.
         */

        typename Trait<TraitType>::Gas *
        findGas(const std::string &name) {
            auto it = std::find_if(gasses.begin(), gasses.end(),
                                [&name](typename Trait<TraitType>::Gas *gas) { return *gas == name; });
            return it == gasses.end() ? nullptr : *it;
        }

        /* -- findState --
         * Tries to find a state in the gas mixture specified by a supplied name. If the
         * state is present, a pointer to this state is returned, otherwise a nullptr is
         * returned.
         *
         * Note that when entry describes a series of states (thus it contains a wildcard
         * character) then a pointer to the first sibling is returned.
         */

        typename Trait<TraitType>::State *
        findState(const StateEntry &entry) {
            // Find gas.

            auto *gas = findGas(entry.gasName);

            if (gas == nullptr)
                return nullptr;

            if (entry.e == "*" && entry.level == electronic) {
                auto &states = entry.charge.empty() ? gas->stateTree : gas->ionicStates;

                if (states.empty())
                    return nullptr;

                return states[0];
            }

            // Find electronic state.

            auto *state = gas->find(entry);

            if (state == nullptr)
                return nullptr;

            if (entry.level == electronic)
                return state;

            if (entry.v == "*" && entry.level == vibrational) {
                if (state->children.empty())
                    return nullptr;

                return state->children[0];
            }

            // Find vibrational state.

            state = state->find(entry);// findState(state, entry);

            if (state == nullptr)
                return nullptr;

            if (entry.level == vibrational)
                return state;

            if (entry.J == "*" && entry.level == rotational) {
                if (state->children.empty())
                    return nullptr;

                return state->children[0];
            }

            // Find rotational state

            state = state->find(entry);//findState(state, entry);

            if (state == nullptr)
                return nullptr;

            if (entry.level == rotational)
                return state;

            return nullptr;
        }

        /* -- addGas --
         * Tries to add a new gas based on a given name and returns a pointer to it. If a
         * gas with the same name already exists it returns a pointer to that gas.
         */

        typename Trait<TraitType>::Gas *
        addGas(const std::string &name) {
            auto *gas = findGas(name);

            if (gas == nullptr)
                return gasses.emplace_back(new typename Trait<TraitType>::Gas(name));

            return gas;

//            auto it = std::find_if(gasses.begin(), gasses.end(), [&name](typename Trait<TraitType>::Gas *gas) {
//                return (*gas == name);
//            });
//
//            if (it == gasses.end()) {
//                return gasses.emplace_back(new typename Trait<TraitType>::Gas(name));
//            }
//
//            return *it;
        }

        /* -- addState --
         * This overload accepts only a reference to a StateEntry object. This is the main
         * function to call when a new state is to be added to the mixture. It will add the
         * corresponding gas and ancestor states if the do not yet exist, before adding the
         * state itself. If the state already exists, then a pointer to this state is
         * returned instead.
         */

        typename Trait<TraitType>::State *
        addState(const StateEntry &entry) {
            auto *gas = addGas(entry.gasName);

            auto *state = addState(gas, entry);

            for (uint8_t lvl = electronic; lvl < entry.level; ++lvl) {
                state = addState(state, entry);
            }

            return state;
        }

        /* -- createCollision --
         * Creates a collision based on a provided CollisionEntry object. The
         * gasses and states involved in the collision are first created and
         * added to the mixture. Then a Collision object is created and its
         * pointer is returned.
         */
        // arguments: smth. like "He + e", "->", "He + e", "Elastic"
        typename Trait<TraitType>::Collision *
        createCollision(const std::string& lhs, const std::string& sep, const std::string& rhs, const std::string& type) {

            std::vector <StateEntry> entry_reactants, entry_products;
            std::vector <uint16_t> entry_stoiCoeff;
            Enumeration::CollisionType entry_type;
            bool entry_isReverse;

            if (!Parse::entriesFromString(lhs, entry_reactants))
                Log<LXCatError>::Error(lhs);
            if (!Parse::entriesFromString(rhs, entry_products, &entry_stoiCoeff))
                Log<LXCatError>::Error(rhs);
            entry_type = Enumeration::getCollisionType(type);
            entry_isReverse = (sep[0] == '<');

            std::vector<typename Trait<TraitType>::State *> reactants, products;
            std::set<typename Trait<TraitType>::Gas *> targetGasses;

            for (auto &stateEntry : entry_reactants) {
                auto *state = reactants.emplace_back(addState(stateEntry));
                targetGasses.insert(state->gas);
            }

            if (targetGasses.size() != 1)
                Log<Message>::Error("Multiple target gasses in a single collision.");

            for (auto &stateEntry : entry_products) {
                products.emplace_back(addState(stateEntry));
            }

            return new typename Trait<TraitType>::Collision(entry_type, reactants, products,
                                                            entry_stoiCoeff, entry_isReverse);
        }

        typename Trait<TraitType>::Collision *
        createCollision(const json_type& rcnf) {

            std::vector <StateEntry> entry_reactants, entry_products;
            std::vector <uint16_t> entry_stoiCoeff;
            Enumeration::CollisionType entry_type;
            bool entry_isReverse;

            if (!Parse::entriesFromJSON(rcnf.at("lhs"), entry_reactants))
                Log<LXCatError>::Error(rcnf.at("lhs").dump(2));
            if (!Parse::entriesFromJSON(rcnf.at("rhs"), entry_products, &entry_stoiCoeff))
                Log<LXCatError>::Error(rcnf.at("rhs").dump(2));
            entry_type = Enumeration::getCollisionType(rcnf.at("type"));

            entry_isReverse = rcnf.at("reverse_also");

            std::vector<typename Trait<TraitType>::State *> reactants, products;
            std::set<typename Trait<TraitType>::Gas *> targetGasses;

            for (auto &stateEntry : entry_reactants) {
                auto *state = reactants.emplace_back(addState(stateEntry));
                targetGasses.insert(state->gas);
            }

            if (targetGasses.size() != 1)
                Log<Message>::Error("Multiple target gasses in a single collision.");

            for (auto &stateEntry : entry_products) {
                products.emplace_back(addState(stateEntry));
            }

            return new typename Trait<TraitType>::Collision(entry_type, reactants, products,
                                                            entry_stoiCoeff, entry_isReverse);
        }

        /* -- loadGasProperties --
         * Loads the gas properties from the database files specified in the main input
         * file. It is declared virtual such that inherited classes can override it to
         * declare any extra properties that they might introduce. They can then call
         * the base class version to load the properties declared in the base class.
         */

        virtual void loadGasProperties(const GasPropertiesSetup &setup) {

            R_GAS_PROPERTY(mass)
            GAS_PROPERTY(harmonicFrequency)
            GAS_PROPERTY(anharmonicFrequency)
            GAS_PROPERTY(electricQuadrupoleMoment)
            GAS_PROPERTY(rotationalConstant)

            // Parse fractions

            const std::regex r(R"(([\w\d]*)\s*=\s*(\d*\.?\d*))");
            std::smatch m;

            for (const auto &fractionStr : setup.fraction) {
                if (!std::regex_search(fractionStr, m, r))
                    Log<Message>::Error("Could not parse gas fractions.");

                const auto &name = m.str(1);

                auto it = std::find_if(gasses.begin(), gasses.end(), [&name](typename Trait<TraitType>::Gas *gas) {
                    return (*gas == name);
                });

                if (it == gasses.end())
                    Log<Message>::Error("Trying to set fraction for non-existent gas: " + name + '.');

                std::stringstream ss(m.str(2));

                if (!(ss >> (*it)->fraction))
                    Log<Message>::Error("Could not parse gas fractions.");

            }

            checkGasFractions();
        }
        virtual void loadGasProperties(const json_type &cnf) {

            R_GAS_PROPERTY_JSON(cnf,mass)
            GAS_PROPERTY_JSON(cnf,harmonicFrequency)
            GAS_PROPERTY_JSON(cnf,anharmonicFrequency)
            GAS_PROPERTY_JSON(cnf,electricQuadrupoleMoment)
            GAS_PROPERTY_JSON(cnf,rotationalConstant)

            // Parse fractions

            const std::regex r(R"(([\w\d]*)\s*=\s*(\d*\.?\d*))");
            std::smatch m;

            for (const std::string &fractionStr : cnf.at("fraction")) {
                if (!std::regex_search(fractionStr, m, r))
                    Log<Message>::Error("Could not parse gas fractions.");

                const auto &name = m.str(1);

                auto it = std::find_if(gasses.begin(), gasses.end(), [&name](typename Trait<TraitType>::Gas *gas) {
                    return (*gas == name);
                });

                if (it == gasses.end())
                    Log<Message>::Error("Trying to set fraction for non-existent gas: " + name + '.');

                std::stringstream ss(m.str(2));

                if (!(ss >> (*it)->fraction))
                    Log<Message>::Error("Could not parse gas fractions.");

            }

            checkGasFractions();
        }

        /* -- checkGasFractions --
         * Checks whether the sum of the gas fractions is equal to 1.
         */

        void checkGasFractions() {
            double norm = 0;

            for (const auto *gas : gasses) {
                norm += gas->fraction;
            }

            if (std::abs(norm - 1.) > 10. * std::numeric_limits<double>::epsilon())
                Log<Message>::Error("Gas fractions are not properly normalized.");
        }

        void checkPopulations() {
            for (auto *gas : gasses)
                gas->checkPopulations();
        }

        /* -- loadStateProperties --
         * Loads the state properties as specified in the input file. It also calls
         * checkPopulations.
         */

        virtual void loadStateProperties(const StatePropertiesSetup &setup,
                                         const WorkingConditions *workingConditions) {

            loadStateProperty(setup.energy, StatePropertyType::energy, workingConditions);
            loadStateProperty(setup.statisticalWeight, StatePropertyType::statisticalWeight, workingConditions);
            loadStateProperty(setup.population, StatePropertyType::population, workingConditions);

            checkPopulations();
        }
        virtual void loadStateProperties(const json_type &cnf,
                                         const WorkingConditions *workingConditions) {

            loadStateProperty(cnf.at("energy"), StatePropertyType::energy, workingConditions);
            loadStateProperty(cnf.at("statisticalWeight"), StatePropertyType::statisticalWeight, workingConditions);
            loadStateProperty(cnf.at("population"), StatePropertyType::population, workingConditions);

            checkPopulations();
        }

        /* -- loadStateProperty --
         * Loads the data concerning a single property of the states (energy, statistical
         * weight, or population) from a vector of entries as supplied by the setup object.
         * First it determines whether the current entry requires loading by direct value,
         * file or function, then it acts accordingly. The property that needs to be set
         * is defined by the propertyType variable. Furthermore, it also needs a pointer
         * to the WorkingConditions structure to access the argument map (which maps its
         * member variables to names by which they are addressed in the input files).
         */

        void loadStateProperty(const std::vector<std::string> &entryVector, StatePropertyType propertyType,
                               const WorkingConditions *workingConditions) {

            for (const auto &line : entryVector) {
                std::string valueString;
                StatePropertyDataType dataType = Parse::statePropertyDataType(line, valueString);

                if (dataType != StatePropertyDataType::file) {
                    StateEntry entry = Parse::propertyStateFromString(line);

                    if (entry.level != none) {
                        auto *state = findState(entry);

                        if (state == nullptr) {
                            Log<PropertyStateError>::Error(entry);
                        }

                        if (dataType == StatePropertyDataType::direct) {
                            double value;

                            if (!Parse::getValue(valueString, value))
                                Log<PropertyValueParseError>::Error(valueString);

                            if (entry.hasWildCard()) {
                                PropertyFunctions::constantValue<TraitType>(state->siblings(), value, propertyType);
                            } else {
                                std::vector<typename Trait<TraitType>::State *> states{state};
                                PropertyFunctions::constantValue<TraitType>(states, value, propertyType);
                            }
                        } else {
                            std::vector<double> arguments;
                            std::string functionName, argumentString;

                            if (!Parse::propertyFunctionAndArguments(valueString, functionName, argumentString))
                                Log<PropertyFunctionParseError>::Error(valueString);

                            if (!Parse::argumentsFromString(argumentString, arguments, workingConditions->argumentMap))
                                Log<PropertyArgumentsError>::Error(argumentString);

                            if (entry.hasWildCard()) {
                                PropertyFunctions::callByName<TraitType>(functionName, state->siblings(),
                                                                         arguments, propertyType);
                            } else {
                                std::vector<typename Trait<TraitType>::State *> states{state};
                                PropertyFunctions::callByName<TraitType>(functionName, states, arguments, propertyType);
                            }
                        }

                    }
                } else {
                    std::vector<std::pair<StateEntry, double>> entries;

                    if (!Parse::statePropertyFile(valueString, entries))
                        Log<FileError>::Error(valueString);

                    for (auto &entry : entries) {
                        auto *state = findState(entry.first);

                        if (state == nullptr)
                            Log<PropertyStateError>::Error(entry.first);

                        if (entry.first.hasWildCard()) {
                            PropertyFunctions::constantValue<TraitType>(state->siblings(), entry.second, propertyType);
                        } else {
                            PropertyFunctions::setStateProperty<TraitType>(state, entry.second, propertyType);
                        }
                    }
                }
            }
        }

        /* -- evaluateStateDensities --
         * Evaluates the densities of all states in the state tree.
         */

        void evaluateStateDensities() {
            for (auto *gas : gasses)
                gas->evaluateStateDensities();
        }

    private:

        /* -- addState --
         * This overload accepts a pointer to a gas. It will check if the electronic
         * ancestor of the state described by the given StateEntry already exists.
         * If it does, a pointer is returned to this state, and if it does not, the
         * ancestor is added its pointer returned.
         */

        typename Trait<TraitType>::State *
        addState(typename Trait<TraitType>::Gas *gas, const StateEntry &entry) {
            auto &states = (entry.charge.empty() ? gas->stateTree : gas->ionicStates);

            auto *state = gas->find(entry);

            if (state == nullptr) {
                state = gas->states.emplace_back(new typename Trait<TraitType>::State(entry, gas));
                return states.emplace_back(state);
            }

            return state;
        }

        /* -- addState --
         * This overload accepts a pointer to a state. It will check if one of the
         * current children of this state is an ancestor of (or equal to) the state
         * described by the given StateEntry. If such a child exists, its pointer is
         * returned, if it does not exist then the appropriate state is created and
         * its pointer returned.
         */

        typename Trait<TraitType>::State *
        addState(typename Trait<TraitType>::State *parent, const StateEntry &entry) {

            auto *state = parent->find(entry);

            if (state == nullptr) {
                state = parent->gas->states.emplace_back(
                        new typename Trait<TraitType>::State(entry, parent->gas, parent));
                return parent->children.emplace_back(state);
            }

            return state;
        }
    };
}


#endif //LOKI_CPP_GASMIXTURE_H
