#include "LoKI-B/StateEntry.h"
#include "LoKI-B/Log.h"
#include <regex>
#include <stdexcept>
#include <string>
#include <lxcat-core.h>

namespace loki
{

StateEntry StateEntry::electronEntry()
{
    static StateEntry el{"e", StateType::charge, "e", "-1", std::string{}, std::string{}, std::string{}};
    return el;
}

StateEntry::StateEntry() : m_id{std::string{}}, m_level(none)
{
}

StateEntry::StateEntry(const std::string &id, StateType level, const std::string &gasName, const std::string &charge,
                       const std::string &e, const std::string &v, const std::string &J)
    : m_id(id), m_level(level), m_gasName(gasName), m_charge(charge), m_e(e), m_v(v), m_J(J)
{
}

/// \todo Update: charge may also be a wildcard (?)
bool StateEntry::hasWildCard() const
{
    switch (m_level)
    {
    case electronic:
        return (m_e == "*");
    case vibrational:
        return (m_v == "*");
    case rotational:
        return (m_J == "*");
    case none:
        return false;
    default:
        break;
    }
    return false;
}

std::ostream &operator<<(std::ostream &os, const StateEntry &entry)
{
    // special handling of the electron. Just write "e".
    if (entry.m_gasName == "e")
    {
        os << entry.m_gasName;
        return os;
    }
    os << entry.m_gasName << '(';
    if (!entry.m_charge.empty())
        os << entry.m_charge << ',';
    os << entry.m_e;
    if (!entry.m_v.empty())
        os << ",v=" << entry.m_v;
    if (!entry.m_J.empty())
        os << ",J=" << entry.m_J;
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
                    StateType stateType = J.empty() == false   ? rotational
                                          : v.empty() == false ? vibrational
                                          : e.empty() == false ? electronic
                                                               : charge;
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

StateEntry entryFromJSON(const std::string &id, const json_type &cnf)
{
    using namespace lxcat_core;

    // TODO: Does disabling this cause problems when loading multiple datasets,
    //       i.e. having multiple electron species around?
    // this is how it is now done for the electron for legacy input
    // if (cnf.at("type").get<std::string>() == "Electron")
    // {
    //     Log<Message>::Warning("Ignoring state attributes for electrons.");
    //     return StateEntry::electronEntry();
    // }

    const auto ser = serialize_species_id(cnf.dump().c_str());

    // We need to split the charge from the gas name, as loki treats them separately.
    std::string gasName(composition(ser));
    gasName = gasName.substr(0, gasName.find("^"));

    const int charge_int = cnf.at("charge").get<int>();
    const std::string charge_str = charge_int ? std::to_string(charge_int) : std::string{};

    // e,v,J are the strings that are passed to the StateEntry constructor.
    std::string e, v, J;

    const auto ele_count = electronic_count(ser);

    if (ele_count > 0)
    {
        if (ele_count != 1)
        {
            throw std::runtime_error("Exactly one electronic state is expected by LoKI-B.");
        }

        e = get_electronic(ser, 0);

        const auto vib_count = vibrational_count(ser);

        if (vib_count > 0)
        {
            if (vib_count == 1)
            {
                v = get_vibrational(ser, 0);

                const auto rot_count = rotational_count(ser);

                if (rot_count > 0)
                {
                    if (rot_count == 1)
                    {
                        J = get_rotational(ser, 0);
                    }
                    else
                    {
                        // For "v" and "J" we assume that the entries form a continuous value-range.
                        std::set<unsigned> J_vals;
                        for (unsigned index = 0; index < rot_count; index++)
                        {
                            J_vals.insert(std::stoi(get_rotational(ser, index)));
                        }
                        if (J_vals.size() != rot_count)
                        {
                            throw std::runtime_error("Duplicate J entries encountered.");
                        }
                        int nJ = *J_vals.rbegin() + 1 - *J_vals.begin();
                        if (nJ != int(J_vals.size()))
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

                if (rotational_count(ser) > 0)
                {
                    throw std::runtime_error("Rotational state identifiers are not allowed when "
                                             "multiple vibrational states are specified.");
                }
                for (unsigned index = 0; index < vib_count; index++)
                {
                    v_vals.insert(std::stoi(get_vibrational(ser, index)));
                }
                // For "v" and "J" we assume that the entries form a continuous value-range.
                if (v_vals.size() != vib_count)
                {
                    throw std::runtime_error("Duplicate v entries encountered.");
                }
                int nv = *v_vals.rbegin() + 1 - *v_vals.begin();
                if (nv != int(v_vals.size()))
                {
                    throw std::runtime_error("Expected a contiguous v-range.");
                }
                v = std::to_string(*v_vals.begin()) + '-' + std::to_string(*v_vals.rbegin());
            }
        }
    }
    StateType stateType = J.empty() == false   ? rotational
                          : v.empty() == false ? vibrational
                          : e.empty() == false ? electronic
                                               : charge;
    return StateEntry{id, stateType, gasName, charge_str, e, v, J};
}

StateEntry propertyStateFromString(const std::string &propertyString)
{
    static const std::regex reState(
        R"(([A-Za-z][A-Za-z0-9]*)\(([-\+]?)\s*,?\s*([-\+'\[\]/\w\*_\^]+)\s*(?:,\s*v\s*=\s*([-\+\w\*]+))?\s*(?:,\s*J\s*=\s*([-\+\d\*]+))?\s*)");
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

} // namespace loki
