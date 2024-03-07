#include "LoKI-B/GasMixture.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/Parse.h"
#include "LoKI-B/PropertyFunctions.h"

#include <filesystem>
#include <regex>

namespace loki
{

GasMixture::~GasMixture()
{
}

Gas *GasMixture::addGas(const std::string& name)
{
    if (findGas(name))
    {
        throw std::logic_error("Attempt to register gas with name '" + name + "' twice.");
    }
    return m_gases.emplace_back(new Gas(name)).get();
}

Gas *GasMixture::ensureGas(const std::string &name)
{
    Gas *gas = findGas(name);
    if (!gas)
    {
        gas = addGas(name);
    }
    return gas;
}

Gas::State *GasMixture::ensureState(const StateEntry &entry)
{
    typename StateMap::iterator it = m_states.find(entry.m_id);
    if (it != m_states.end())
    {
        return it->second;
    }
    Gas::State *state = ensureGas(entry.m_gasName)->ensureState(entry);
    m_states[entry.m_id] = state;
    return state;
}

Gas::State *GasMixture::findStateById(const std::string &stateId)
{
    typename StateMap::iterator it = m_states.find(stateId);
    return it == m_states.end() ? nullptr : it->second;
}
void GasMixture::print(std::ostream &os)
{
    for (const auto &gas : m_gases)
    {
        os << "Gas: " << gas->name() << std::endl;
        gas->print(os);
    }
}

void GasMixture::checkGasFractions()
{
    double norm = 0;

    for (const auto &gas : m_gases)
    {
        norm += gas->fraction;
    }

    if (std::abs(norm - 1.) > 10. * std::numeric_limits<double>::epsilon())
        Log<Message>::Error("Gas fractions are not properly normalized.");
}

void GasMixture::checkPopulations()
{
    for (auto &gas : m_gases)
        gas->checkPopulations();
}

Gas *GasMixture::findGas(const std::string &name)
{
    auto it = std::find_if(m_gases.begin(), m_gases.end(),
                           [&name](const std::unique_ptr<Gas> &gas) { return gas->name() == name; });
    return it == m_gases.end() ? nullptr : it->get();
}

Gas::State *GasMixture::findState(const StateEntry &entry)
{
    Gas *gas = findGas(entry.m_gasName);
    return gas ? gas->findState(entry) : nullptr;
}

Gas::State::ChildContainer GasMixture::findStates(const StateEntry &entry)
{
    using ChildContainer = Gas::State::ChildContainer;
    Gas::State* state = findState(entry);
    return (!state) ? ChildContainer{} : entry.hasWildCard() ? state->siblings() : ChildContainer{ state };
}

void GasMixture::loadStateProperty(const std::filesystem::path &basePath, const std::vector<std::string> &entryVector, 
                                   StatePropertyType propertyType, const WorkingConditions *workingConditions)
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
            /** \todo Should this be part of propertyStateFromString?
             *        Is there a reason to accept a 'none'-result?
             */
            if (entry.m_level == none)
            {
                throw std::runtime_error("loadStateProperty: illegal "
                                "state identifier '" + line + "'.");
            }
            Gas::State::ChildContainer states = findStates(entry);
            if (states.empty())
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
                PropertyFunctions::constantValue(states, value, propertyType);
            }
            else
            {
                // expr is not a number. We treat is a function (maybe with arguments).
                static const std::regex reFuncArgs(R"(\s*(\w+)@?(.*))");
                std::smatch fm;
                if (!std::regex_match(expr, fm, reFuncArgs))
                {
                    throw std::runtime_error("Could not parse function "
                        "name and argument list from string '"
                        + expr + "'.");
                }
                const std::string functionName = fm.str(1);
                const std::string argumentString = fm.str(2);

                // create an argument list for the function (possibly empty)
                std::vector<double> arguments;
                static const std::regex reArgList(R"(\s*([\w\.]+)\s*(?:[,\]]|$))");
                for (auto it = std::sregex_iterator(argumentString.begin(), argumentString.end(), reArgList);
                     it != std::sregex_iterator(); ++it)
                {
                    const std::string arg{it->str(1)};
                    double pvalue;
                    if (Parse::getValue(arg, pvalue))
                    {
                        // pvalue set by getValue
                    }
                    else if (workingConditions->findParameter(arg))
                    {
                        pvalue = workingConditions->getParameter(arg);
                    }
                    else
                    {
                        throw std::runtime_error("Argument '" + arg + "' is neither a numerical "
                                    "pvalue, nor a known parameter name.");
                    }
                    arguments.emplace_back(pvalue);
                }
                PropertyFunctions::callByName(functionName, states, arguments, propertyType);
            }
        }
        else
        {
            std::vector<std::pair<StateEntry, double>> entries;
            // const std::string fileName = line;
            std::filesystem::path fileName(line);

            if (fileName.is_relative()) {
                fileName = basePath.parent_path() / fileName;
            }

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

            /** \todo Also in case we read a property file, the expression type could be a function.
             *        It would be nice if the external file case behaves the same as the inline case.
             *        This could be achieved by first assembling a collection of tasks (either from
             *        the settings node or from file), then execute these tasks. That will be more
             *        general and simplify this function at the same time.
             */
            for (auto &entry : entries)
            {
                Gas::State::ChildContainer states = findStates(entry.first);
                if (states.empty())
                {
                    throw std::runtime_error("loadStateProperty: could not find "
                                    "state or state group '" + entry.first.m_id + "'.");
                }
                PropertyFunctions::constantValue(states, entry.second, propertyType);
            }
        }
    }
}

void GasMixture::evaluateReducedDensities()
{
    for (auto &gas : m_gases)
    {
        gas->evaluateReducedDensities();
    }
}

void GasMixture::loadStateProperties(const std::filesystem::path &basePath, const StatePropertiesSetup &setup, const WorkingConditions *workingConditions)
{

    loadStateProperty(basePath, setup.energy, StatePropertyType::energy, workingConditions);
    loadStateProperty(basePath, setup.statisticalWeight, StatePropertyType::statisticalWeight, workingConditions);
    loadStateProperty(basePath, setup.population, StatePropertyType::population, workingConditions);

    checkPopulations();
}

void GasMixture::loadStateProperties(const std::filesystem::path &basePath, const json_type &cnf, const WorkingConditions *workingConditions)
{

    loadStateProperty(basePath, cnf.at("energy"), StatePropertyType::energy, workingConditions);
    loadStateProperty(basePath, cnf.at("statisticalWeight"), StatePropertyType::statisticalWeight, workingConditions);
    loadStateProperty(basePath, cnf.at("population"), StatePropertyType::population, workingConditions);

    checkPopulations();
}

void GasMixture::loadGasProperties(const std::filesystem::path &basePath, const GasPropertiesSetup &setup)
{
    readGasPropertyFile(basePath, m_gases, setup.mass, "mass", true,
        [](Gas& gas, double value) { gas.mass=value; } );
    readGasPropertyFile(basePath, m_gases, setup.harmonicFrequency, "harmonicFrequency", false,
        [](Gas& gas, double value) { gas.harmonicFrequency=value; } );
    readGasPropertyFile(basePath, m_gases, setup.anharmonicFrequency, "anharmonicFrequency", false,
        [](Gas& gas, double value) { gas.anharmonicFrequency=value; } );
    readGasPropertyFile(basePath, m_gases, setup.electricQuadrupoleMoment, "electricQuadrupoleMoment", false,
        [](Gas& gas, double value) { gas.electricQuadrupoleMoment=value; } );
    readGasPropertyFile(basePath, m_gases, setup.rotationalConstant, "rotationalConstant", false,
        [](Gas& gas, double value) { gas.rotationalConstant=value; } );

    // Parse fractions
    const std::regex r(R"(([\w\d]*)\s*=\s*(\d*\.?\d*))");
    std::smatch m;

    for (const auto &fractionStr : setup.fraction)
    {
        if (!std::regex_search(fractionStr, m, r))
            Log<Message>::Error("Could not parse gas fractions.");

        const auto &name = m.str(1);

        auto it = std::find_if(m_gases.begin(), m_gases.end(),
                               [&name](std::unique_ptr<Gas> &gas) { return (gas->name() == name); });

        if (it == m_gases.end())
            Log<Message>::Error("Trying to set fraction for non-existent gas: " + name + '.');

        std::stringstream ss(m.str(2));

        if (!(ss >> (*it)->fraction))
            Log<Message>::Error("Could not parse gas fractions.");
    }

    checkGasFractions();
}

void GasMixture::loadGasProperties(const std::filesystem::path &basePath, const json_type &cnf)
{
    readGasProperty(basePath, m_gases, cnf, "mass", true,
        [](Gas& gas, double value) { gas.mass=value; } );
    readGasProperty(basePath, m_gases, cnf, "harmonicFrequency", false,
        [](Gas& gas, double value) { gas.harmonicFrequency=value; } );
    readGasProperty(basePath, m_gases, cnf, "anharmonicFrequency", false,
        [](Gas& gas, double value) { gas.anharmonicFrequency=value; } );
    readGasProperty(basePath, m_gases, cnf, "electricQuadrupoleMoment", false,
        [](Gas& gas, double value) { gas.electricQuadrupoleMoment=value; } );
    readGasProperty(basePath, m_gases, cnf, "rotationalConstant", false,
        [](Gas& gas, double value) { gas.rotationalConstant=value; } );

    // Parse fractions
    const std::regex r(R"(([\w\d]*)\s*=\s*(\d*\.?\d*))");
    std::smatch m;

    for (const std::string fractionStr : cnf.at("fraction"))
    {
        if (!std::regex_search(fractionStr, m, r))
            Log<Message>::Error("Could not parse gas fractions.");
        const auto &name = m.str(1);
        auto it = std::find_if(m_gases.begin(), m_gases.end(),
                               [&name](std::unique_ptr<Gas> &gas) { return (gas->name() == name); });
        if (it == m_gases.end())
            Log<Message>::Error("Trying to set fraction for non-existent gas: " + name + '.');
        std::stringstream ss(m.str(2));

        if (!(ss >> (*it)->fraction))
            Log<Message>::Error("Could not parse gas fractions.");
    }
    checkGasFractions();
}

} // namespace loki
