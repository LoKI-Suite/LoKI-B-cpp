#include "LoKI-B/GasMixtureBase.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/Parse.h"
#include "LoKI-B/PropertyFunctions.h"
/// \todo Remove me once INPUT is no longer used.
#include "LoKI-B/StandardPaths.h"


#include <regex>

namespace loki
{

GasBase *GasMixtureBase::addGas(GasBase *gas)
{
    if (findGas(gas->name))
    {
        throw std::logic_error("Attempt to register gas with name '" + gas->name + "' twice.");
    }
    return m_gases.emplace_back(gas).get();
}

void GasMixtureBase::print(std::ostream &os)
{
    for (const auto &gas : m_gases)
    {
        os << "Gas: " << gas->name << std::endl;
        gas->print(os);
    }
}

void GasMixtureBase::checkGasFractions()
{
    double norm = 0;

    for (const auto &gas : m_gases)
    {
        norm += gas->fraction;
    }

    if (std::abs(norm - 1.) > 10. * std::numeric_limits<double>::epsilon())
        Log<Message>::Error("Gas fractions are not properly normalized.");
}

void GasMixtureBase::checkPopulations()
{
    for (auto &gas : m_gases)
        gas->checkPopulations();
}

GasBase *GasMixtureBase::findGas(const std::string &name)
{
    auto it = std::find_if(m_gases.begin(), m_gases.end(),
                           [&name](const std::unique_ptr<GasBase> &gas) { return gas->name == name; });
    return it == m_gases.end() ? nullptr : it->get();
}

GasBase::StateBase *GasMixtureBase::findState(const StateEntry &entry)
{
    GasBase *gas = findGas(entry.gasName);
    return gas ? gas->findState(entry) : nullptr;
}

void GasMixtureBase::loadStateProperty(const std::vector<std::string> &entryVector, StatePropertyType propertyType,
                                       const WorkingConditions *workingConditions)
{

    for (const auto &line : entryVector)
    {
        // look for a line of the form "S = E"
        static const std::regex reProperty(R"((\S+)\s+=\s+(\S+)\s*$)");
        std::smatch m;
        if (std::regex_search(line, m, reProperty))
        {
            // Found. E can be a literal value, or a function with arguments.
            // value or function
            const std::string state_id = m.str(1);
            const std::string expr = m.str(2);

            // 1. get the state or state group that the expression will
            //    be applied to.
            const StateEntry entry = propertyStateFromString(state_id);
            if (entry.level == none)
            {
                throw std::runtime_error("loadStateProperty: illegal "
                                "state identifier '" + line + "'.");
            }
            GasBase::StateBase *state = findState(entry);
            if (state == nullptr)
            {
                throw std::runtime_error("loadStateProperty: could not find "
                                "state or state group '" + line + "'.");
            }

            // 2. Now apply the expression.
            // Try to parse expr as a number first...
            double value;
            if (Parse::getValue(expr, value))
            {
                // expr is a number, now parsed into value.
                if (entry.hasWildCard())
                {
                    PropertyFunctions::constantValue(state->siblings(), value, propertyType);
                }
                else
                {
                    /// \todo Can we avoid creation of the intermediate vector?
                    std::vector<GasBase::StateBase *> states{state};
                    PropertyFunctions::constantValue(states, value, propertyType);
                }
            }
            else
            {
                // expr is not a number. We treat is a function (maybe with arguments).
                static const std::regex reFuncArgs(R"(\s*(\w+)@?(.*))");
                std::smatch m;
                if (!std::regex_match(expr, m, reFuncArgs))
                {
                    throw std::runtime_error("Could not parse function "
                        "name and argument list from string '"
                        + expr + "'.");
                }
                const std::string functionName = m.str(1);
                const std::string argumentString = m.str(2);

                // create an argument list for the function (possibly empty)
                std::vector<double> arguments;
                static const std::regex reArgList(R"(\s*([\w\.]+)\s*(?:[,\]]|$))");
                for (auto it = std::sregex_iterator(argumentString.begin(), argumentString.end(), reArgList);
                     it != std::sregex_iterator(); ++it)
                {
                    const std::string arg{it->str(1)};
                    double value;
                    if (Parse::getValue(arg, value))
                    {
                        // value set by getValue
                    }
                    else if (workingConditions->argumentMap().count(arg))
                    {
                        value = *workingConditions->argumentMap().at(arg);
                    }
                    else
                    {
                        throw std::runtime_error("Argument '" + arg + "' is neither a numerical "
                                    "value, nor a known parameter name.");
                    }
                    arguments.emplace_back(value);
                }

                if (entry.hasWildCard())
                {
                    PropertyFunctions::callByName(functionName, state->siblings(), arguments, propertyType);
                }
                else
                {
                    /// \todo Can we avoid creation of the intermediate vector?
                    std::vector<GasBase::StateBase *> states{state};
                    PropertyFunctions::callByName(functionName, states, arguments, propertyType);
                }
            }
        }
        else
        {
            std::vector<std::pair<StateEntry, double>> entries;
            const std::string fileName = INPUT "/" + line;

            try
            {
                statePropertyFile(fileName, entries);
            }
            catch (std::exception& exc)
            {
                throw std::runtime_error("Error configuring state property '"
                                        + statePropertyName(propertyType) + "': "
                                        + std::string{exc.what()} );
            }

            for (auto &entry : entries)
            {
                GasBase::StateBase *state = findState(entry.first);

                if (state == nullptr)
                    Log<PropertyStateError>::Error(entry.first);

                if (entry.first.hasWildCard())
                {
                    PropertyFunctions::constantValue(state->siblings(), entry.second, propertyType);
                }
                else
                {
                    PropertyFunctions::setStateProperty(state, entry.second, propertyType);
                }
            }
        }
    }
}

void GasMixtureBase::evaluateStateDensities()
{
    for (auto &gas : m_gases)
        gas->evaluateStateDensities();
}

void GasMixtureBase::loadStateProperties(const StatePropertiesSetup &setup, const WorkingConditions *workingConditions)
{

    loadStateProperty(setup.energy, StatePropertyType::energy, workingConditions);
    loadStateProperty(setup.statisticalWeight, StatePropertyType::statisticalWeight, workingConditions);
    loadStateProperty(setup.population, StatePropertyType::population, workingConditions);

    checkPopulations();
}

void GasMixtureBase::loadStateProperties(const json_type &cnf, const WorkingConditions *workingConditions)
{

    loadStateProperty(cnf.at("energy"), StatePropertyType::energy, workingConditions);
    loadStateProperty(cnf.at("statisticalWeight"), StatePropertyType::statisticalWeight, workingConditions);
    loadStateProperty(cnf.at("population"), StatePropertyType::population, workingConditions);

    checkPopulations();
}

void GasMixtureBase::loadGasProperties(const GasPropertiesSetup &setup)
{
    readGasPropertyFile(m_gases, setup.mass, "mass", true,
        [](GasBase& gas, double value) { gas.mass=value; } );
    readGasPropertyFile(m_gases, setup.harmonicFrequency, "harmonicFrequency", false,
        [](GasBase& gas, double value) { gas.harmonicFrequency=value; } );
    readGasPropertyFile(m_gases, setup.anharmonicFrequency, "anharmonicFrequency", false,
        [](GasBase& gas, double value) { gas.anharmonicFrequency=value; } );
    readGasPropertyFile(m_gases, setup.electricQuadrupoleMoment, "electricQuadrupoleMoment", false,
        [](GasBase& gas, double value) { gas.electricQuadrupoleMoment=value; } );
    readGasPropertyFile(m_gases, setup.rotationalConstant, "rotationalConstant", false,
        [](GasBase& gas, double value) { gas.rotationalConstant=value; } );

    // Parse fractions
    const std::regex r(R"(([\w\d]*)\s*=\s*(\d*\.?\d*))");
    std::smatch m;

    for (const auto &fractionStr : setup.fraction)
    {
        if (!std::regex_search(fractionStr, m, r))
            Log<Message>::Error("Could not parse gas fractions.");

        const auto &name = m.str(1);

        auto it = std::find_if(m_gases.begin(), m_gases.end(),
                               [&name](std::unique_ptr<GasBase> &gas) { return (gas->name == name); });

        if (it == m_gases.end())
            Log<Message>::Error("Trying to set fraction for non-existent gas: " + name + '.');

        std::stringstream ss(m.str(2));

        if (!(ss >> (*it)->fraction))
            Log<Message>::Error("Could not parse gas fractions.");
    }

    checkGasFractions();
}

void GasMixtureBase::loadGasProperties(const json_type &cnf)
{
    readGasProperty(m_gases, cnf, "mass", true,
        [](GasBase& gas, double value) { gas.mass=value; } );
    readGasProperty(m_gases, cnf, "harmonicFrequency", false,
        [](GasBase& gas, double value) { gas.harmonicFrequency=value; } );
    readGasProperty(m_gases, cnf, "anharmonicFrequency", false,
        [](GasBase& gas, double value) { gas.anharmonicFrequency=value; } );
    readGasProperty(m_gases, cnf, "electricQuadrupoleMoment", false,
        [](GasBase& gas, double value) { gas.electricQuadrupoleMoment=value; } );
    readGasProperty(m_gases, cnf, "rotationalConstant", false,
        [](GasBase& gas, double value) { gas.rotationalConstant=value; } );

    // Parse fractions
    const std::regex r(R"(([\w\d]*)\s*=\s*(\d*\.?\d*))");
    std::smatch m;

    for (const std::string &fractionStr : cnf.at("fraction"))
    {
        if (!std::regex_search(fractionStr, m, r))
            Log<Message>::Error("Could not parse gas fractions.");
        const auto &name = m.str(1);
        auto it = std::find_if(m_gases.begin(), m_gases.end(),
                               [&name](std::unique_ptr<GasBase> &gas) { return (gas->name == name); });
        if (it == m_gases.end())
            Log<Message>::Error("Trying to set fraction for non-existent gas: " + name + '.');
        std::stringstream ss(m.str(2));

        if (!(ss >> (*it)->fraction))
            Log<Message>::Error("Could not parse gas fractions.");
    }
    checkGasFractions();
}

} // namespace loki
