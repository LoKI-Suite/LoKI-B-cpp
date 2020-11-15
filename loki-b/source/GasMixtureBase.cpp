#include "LoKI-B/GasMixtureBase.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/Parse.h"
#include "LoKI-B/PropertyFunctions.h"

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
        std::string valueString;
        StatePropertyDataType dataType = Parse::statePropertyDataType(line, valueString);

        if (dataType != StatePropertyDataType::file)
        {
            StateEntry entry = propertyStateFromString(line);

            if (entry.level != none)
            {
                GasBase::StateBase *state = findState(entry);

                if (state == nullptr)
                {
                    Log<PropertyStateError>::Error(entry);
                }

                if (dataType == StatePropertyDataType::direct)
                {
                    double value;

                    if (!Parse::getValue(valueString, value))
                        Log<PropertyValueParseError>::Error(valueString);

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
                    std::vector<double> arguments;
                    std::string functionName, argumentString;

                    if (!Parse::propertyFunctionAndArguments(valueString, functionName, argumentString))
                        Log<PropertyFunctionParseError>::Error(valueString);

                    if (!Parse::argumentsFromString(argumentString, arguments, workingConditions->argumentMap))
                        Log<PropertyArgumentsError>::Error(argumentString);

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
        }
        else
        {
            std::vector<std::pair<StateEntry, double>> entries;
            const std::string fileName = INPUT "/" + valueString;

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
