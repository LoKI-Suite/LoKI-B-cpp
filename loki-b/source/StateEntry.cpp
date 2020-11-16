#include "LoKI-B/StateEntry.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/Parse.h"
#include "LoKI-B/StandardPaths.h"
#include <fstream>
#include <stdexcept>
#include <regex>

namespace loki
{

StateEntry StateEntry::electronEntry()
{
    static StateEntry el{"e", StateType::charge, "e", "-", std::string{}, std::string{}, std::string{}};
    return el;
}

StateEntry::StateEntry() : m_id{std::string{}}, level(none)
{
}

StateEntry::StateEntry(const std::string &id, StateType level, const std::string &gasName, const std::string &charge,
                       const std::string &e, const std::string &v, const std::string &J)
    : m_id(id), level(level), charge(charge), gasName(gasName), e(e), v(v), J(J)
{
}

/// \todo Update: charge may also be a wildcard (?)
bool StateEntry::hasWildCard() const
{
    switch (level)
    {
    case electronic:
        return (e == "*");
    case vibrational:
        return (v == "*");
    case rotational:
        return (J == "*");
    case none:
        return false;
    }
    return false;
}

std::ostream &operator<<(std::ostream &os, const StateEntry &entry)
{
    // special handling of the electron. Just write "e".
    if (entry.gasName == "e")
    {
        os << entry.gasName;
        return os;
    }
    os << entry.gasName << '(';
    if (!entry.charge.empty())
        os << entry.charge << ',';
    os << entry.e;
    if (!entry.v.empty())
        os << ",v=" << entry.v;
    if (!entry.J.empty())
        os << ",J=" << entry.J;
    os << ')';

    return os;
}

namespace
{

// Parenthesize a string and return the result. This makes a regex group
// and is used for readability.
std::string group(const std::string &term)
{
    return "(" + term + ")";
}

} // namespace

void entriesFromString(const std::string stateString, std::vector<StateEntry> &entries,
                       std::vector<uint16_t> *stoiCoeff)
{
    const std::string plus_sign = "\\+";
    const std::string ws = "\\s*";
    const std::string coef = group("\\d*");
    const std::string state_id = group(".*?");
    const std::string gasname = group("[A-Z][[:alnum:]]*");
    const std::string electron = group("e");
    const std::string gas_particle = gasname + "\\(" + state_id + "\\)";
    const std::string particle = electron + "|" + gas_particle;
    const std::string term = ws + plus_sign + ws + coef + ws + group(particle) + ws;

    // we expect smth. like 'term + term + term'. Change the target string into
    // '+ term + term + term', so we can parse plus-term pairs and need no special
    // measures for the missing leading '+'
    const std::string parseString = "+ " + stateString;
    std::string remainder = parseString;

    const std::regex term_expr('^' + term);
    std::smatch res;
    while (remainder.size() && std::regex_search(remainder, res, term_expr))
    {
        const std::string part_remainder = remainder;
        try
        {
            // as a result of the special handling of the electron, two types of
            // submatch-sequences are possible, as shown below. Also the submatch
            // indices are indicated.
            //       0        1       2         3   4      5
            // a. <group> [  <coef> <group> [  'e'  <>     <>    ] ]
            // b. <group> [  <coef> <group> [  <>  <gas> <state> ] ]

            // the stoichiometric coefficient (empty corresponds to 1).
            const std::string c = res[1];
            const std::string id = res.str(2);
            std::string g, s;
            // state id components:
            std::string q, e, v, J;
            if (res.str(3).empty())
            {
                g = res[4];
                s = res[5];

                try
                {

                    // now parse the state id list
                    std::smatch state_res;
                    std::string stateRemainder = s;
                    const std::regex charge_expr{"^" + group("\\+*|\\-*") + ","};
                    if (std::regex_search(stateRemainder, state_res, charge_expr))
                    {
                        q = state_res[1];
                        stateRemainder = state_res.suffix().str();
                    }
                    const std::regex elec_expr{"^" + group("[^,=]+") + ",?"};
                    if (std::regex_search(stateRemainder, state_res, elec_expr))
                    {
                        e = state_res[1];
                        stateRemainder = state_res.suffix().str();
                    }
                    else
                    {
                        throw std::runtime_error("Bad electronic id, starting at '" + stateRemainder + ".");
                    }
                    const std::regex vib_expr{"^v=" + group("[^,=]+") + ",?"};
                    if (!stateRemainder.empty())
                    {
                        if (std::regex_search(stateRemainder, state_res, vib_expr))
                        {
                            v = state_res[1];
                            stateRemainder = state_res.suffix().str();
                        }
                        else
                        {
                            throw std::runtime_error("Bad vibrational id, starting at '" + stateRemainder + ".");
                        }
                    }
                    const std::regex rot_expr{"^J=" + group("[^,=]+") + ",?"};
                    if (!stateRemainder.empty())
                    {
                        if (std::regex_search(stateRemainder, state_res, rot_expr))
                        {
                            J = state_res[1];
                            stateRemainder = state_res.suffix().str();
                        }
                        else
                        {
                            throw std::runtime_error("Bad rotational id, starting at '" + stateRemainder + ".");
                        }
                    }
                    if (!stateRemainder.empty())
                    {
                        throw std::runtime_error("Found trailing characters '" + stateRemainder + "'.");
                    }
                    StateType stateType =
                        J.empty() == false
                            ? rotational
                            : v.empty() == false ? vibrational : e.empty() == false ? electronic : charge;
                    entries.push_back(StateEntry(id, stateType, g, q, e, v, J));
                    if (stoiCoeff)
                    {
                        stoiCoeff->push_back(c.empty() ? 1 : std::stoi(c));
                    }
                }
                catch (std::exception &exc)
                {
                    throw std::runtime_error("While parsing state identifier list '" + s + "':\n" +
                                             std::string(exc.what()));
                }
            }
            else
            {
                g = res[3];
                s = std::string{};
                if (g != "e")
                {
                    throw std::logic_error("Expected 'e', found '" + g + ".");
                }
                // special handling for the electron:
                entries.push_back(StateEntry::electronEntry());
                if (stoiCoeff)
                {
                    stoiCoeff->push_back(c.empty() ? 1 : std::stoi(c));
                }
            }
#if 0
        std::cout << "stoich = '" << c
            << "', particle = '" << g << "'."
            << "', state id = '" << s
                << "', q = '" << q
                << "', e = '" << e
                << "', v = '" << v
                << "', J = '" << J
                << "'" << std::endl;
#endif
            remainder = res.suffix().str();
        }
        catch (std::exception &exc)
        {
            throw std::runtime_error("While parsing particle '" + part_remainder + "':\n" + std::string{exc.what()});
        }
    }
    // incomplete parse?
    if (!remainder.empty())
    {
        // if parsing fails at the beginning, show the original string without the
        // artificial '+ ' prepended.
        throw std::runtime_error("Parsing of stoichiometric array failed at '" +
                                 (remainder == parseString ? stateString : remainder) + "'.");
    }
}

StateEntry entryFromJSON(const json_type &cnf)
{
    const std::string gasName = cnf.at("particle");
    const std::string id = cnf.at("id");
    // this is how it is now done for the electron for legacy input
    if (gasName == "e")
    {
        Log<Message>::Warning("Ignoring state attributes for electrons.");
        return StateEntry::electronEntry();
    }
    const int charge_int = cnf.at("charge").get<int>();
    const std::string charge_str = charge_int ? std::to_string(charge_int) : std::string{};
    const json_type &el_cnf = cnf.at("electronic");
    if (el_cnf.size() != 1)
    {
        throw std::runtime_error("Exactly one electronic state is expected by LoKI-B.");
    }
    // e,v,J are the strings that are passed to the StateEntry constructor.
    const std::string e = el_cnf[0].at("e");
    std::string v, J;
    if (el_cnf[0].contains("vibrational"))
    {
        const json_type &vib_cnf = el_cnf[0].at("vibrational");
        if (vib_cnf.size() == 0)
        {
            throw std::runtime_error("At least one vibrational state is expected by LoKI-B.");
        }
        else if (vib_cnf.size() == 1)
        {
            // we expect a number, but sometimes a string is encountered, like "10+"
            v = vib_cnf[0].at("v").type() == json_type::value_t::string ? vib_cnf[0].at("v").get<std::string>()
                                                                        : std::to_string(vib_cnf[0].at("v").get<int>());
            if (vib_cnf[0].contains("rotational"))
            {
                const json_type &rot_cnf = vib_cnf[0].at("rotational");
                if (rot_cnf.size() == 0)
                {
                    throw std::runtime_error("At least one rotational state is expected by LoKI-B.");
                }
                else if (rot_cnf.size() == 1)
                {
                    J = std::to_string(rot_cnf[0].at("J").get<int>());
                }
                else
                {
                    // For "v" and "J" we assume that the entries form a continuous value-range.
                    std::set<unsigned> J_vals;
                    for (const auto &Jentry : rot_cnf)
                    {
                        J_vals.insert(Jentry.at("J").get<int>());
                    }
                    if (J_vals.size() != rot_cnf.size())
                    {
                        throw std::runtime_error("Duplicate J entries encountered.");
                    }
                    int nJ = *J_vals.rbegin() + 1 - *J_vals.begin();
                    if (nJ != J_vals.size())
                    {
                        throw std::runtime_error("Expected a contiguous J-range.");
                    }
                    J = std::to_string(*J_vals.begin()) + '-' + std::to_string(*J_vals.rbegin());
                }
            }
        }
        else
        {
            std::set<unsigned> v_vals;
            for (const auto &ventry : vib_cnf)
            {
                if (ventry.contains("rotational"))
                {
                    throw std::runtime_error("Rotational states identifiers are not allowed when "
                                             "multiple virbational states are specified.");
                }
                v_vals.insert(ventry.at("v").get<int>());
            }
            // For "v" and "J" we assume that the entries form a continuous value-range.
            if (v_vals.size() != vib_cnf.size())
            {
                throw std::runtime_error("Duplicate v entries encountered.");
            }
            int nv = *v_vals.rbegin() + 1 - *v_vals.begin();
            if (nv != v_vals.size())
            {
                throw std::runtime_error("Expected a contiguous v-range.");
            }
            v = std::to_string(*v_vals.begin()) + '-' + std::to_string(*v_vals.rbegin());
        }
    }
    StateType stateType =
        J.empty() == false ? rotational : v.empty() == false ? vibrational : e.empty() == false ? electronic : charge;
    return StateEntry{id, stateType, gasName, charge_str, e, v, J};
}

StateEntry propertyStateFromString(const std::string &propertyString)
{
    static const std::regex reState(
        R"(([A-Za-z][A-Za-z0-9]*)\(([-\+]?)\s*,?\s*([-\+'\[\]/\w\*]+)\s*(?:,\s*v\s*=\s*([-\+\w\*]+))?\s*(?:,\s*J\s*=\s*([-\+\d\*]+))?\s*)");
    std::smatch m;

    if (!std::regex_search(propertyString, m, reState))
        return {};

    if (m.str(1).empty() || m.str(3).empty())
        return {};

    StateType stateType;

    if (m.str(4).empty())
    {
        stateType = electronic;
    }
    else if (m.str(5).empty())
    {
        stateType = vibrational;
    }
    else
    {
        stateType = rotational;
    }
    return {m.str(0), stateType, m.str(1), m.str(2), m.str(3), m.str(4), m.str(5)};
}

void statePropertyFile(const std::string &fileName, std::vector<std::pair<StateEntry, double>> &entries)
{
    /** \bug 'S 1.2.3' will be accepted by this regex. Subsequently,
     *        getValue will result in the value 1.2, since that does
     *        not care about trailing characters.
     */
    static const std::regex expr(R"((.*?)\s+([\d\.e+-]+)\s*(?:\n|$))");

    std::string fileBuffer;
    if (!Parse::stringBufferFromFile(fileName, fileBuffer))
    {
        Log<Message>::Error("Could not open state property file '"
                            + fileName + "' for reading.");
    }
    std::stringstream ss{fileBuffer};
    std::string line;
    while (std::getline(ss, line))
    {
        if (line.size() < 3)
        {
            continue;
        }
        std::smatch m;
        if (!std::regex_search(line, m, expr))
        {
            throw std::runtime_error("Syntax error in file '"
                            + fileName + "', line '"
                            + line + "': expected a state name and a numeric argument.");
        }
        const std::string stateString = m.str(1);
        const std::string valueString = m.str(2);
        double value;
        if (!Parse::getValue(valueString, value))
        {
            throw std::runtime_error("Syntax error in file '"
                            + fileName + "', line '"
                            + line + "': could not convert argument '"
                            + valueString + "' to a number.");
        }
        entries.emplace_back(propertyStateFromString(stateString), value);
    }
}

} // namespace loki
