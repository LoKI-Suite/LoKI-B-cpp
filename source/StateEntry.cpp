#include "LoKI-B/StateEntry.h"
#include "LoKI-B/Enumeration.h"
#include "LoKI-B/Log.h"
#include <regex>
#include <set>
#include <stdexcept>

namespace loki
{

StateEntry StateEntry::electronEntry()
{
    static StateEntry el{"e", StateType::charge, "e", "-", std::string{}, std::string{}, std::string{}};
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
    os << entry.m_gasName;

    // special handling of the electron. Just write "e".
    if (entry.m_gasName == "e" || entry.m_level == StateType::root)
    {
        return os;
    }

    os << '(';

    if (!entry.m_charge.empty())
        os << entry.m_charge;

    if (!entry.m_e.empty())
    {
        if (!entry.m_charge.empty())
            os << ',';
        os << entry.m_e;
    }

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

std::string serialize_charge(int charge) {
    if (charge == 0) {
        return "";
    }

    return std::string(std::abs(charge), charge > 0 ? '+' : '-');
}

StateEntry entryFromJSON(const std::string &id, const json_type &cnf)
{
    const json_type &ser_cnf = cnf.at("serialized");
    const auto species_type = cnf.at("detailed").at("type").get<std::string_view>();

    std::string gas_name = ser_cnf.at("composition").at("summary");

    // Split the charge from the gas name.
    const auto pos = gas_name.find("^");
    if (pos != std::string::npos)
    {
        gas_name = gas_name.substr(0, pos);
    }

    const int charge_int = cnf.at("detailed").at("charge").get<int>();
    const std::string charge_str = serialize_charge(charge_int);

    // e,v,J are the strings that are passed to the StateEntry constructor.
    std::string e, v, J;
    if (ser_cnf.contains("electronic"))
    {
        const json_type &el_cnf = ser_cnf.at("electronic");
        if (el_cnf.is_array())
        {
            for (nlohmann::json::size_type i = 0; i < el_cnf.size(); i++) {
                e.append(el_cnf[i].at("summary").get<std::string>());
                if (i < el_cnf.size() - 1) e.append("|");
            }
        } else if (el_cnf.is_object())
        {
            e = el_cnf.at("summary");
            // throw std::runtime_error("Exactly one electronic state is expected by LoKI-B.");
        }
        if (el_cnf.contains("vibrational"))
        {
            const json_type &vib_cnf = el_cnf.at("vibrational");
            if (vib_cnf.is_object())
            {
                v = vib_cnf.at("summary").get<std::string>();

                if (vib_cnf.contains("rotational"))
                {
                    const json_type &rot_cnf = vib_cnf.at("rotational");
                    if (rot_cnf.is_object())
                    {
                        J = rot_cnf.at("summary").get<std::string>();
                    }
                    else
                    {
                        if (rot_cnf.size() == 0)
                        {
                            throw std::runtime_error("At least one rotational state is expected by LoKI-B.");
                        }

                        // For "v" and "J" we assume that the entries form a continuous value-range.
                        std::set<unsigned> J_vals;
                        for (const auto &Jentry : rot_cnf)
                        {
                            // Avoid unexpected behavior for rotational states of molecules with multiple rotational
                            // quanta (e.g. water).
                            const auto &J_str = Jentry.at("summary").get<std::string>();

                            if (!(species_type == "HomonuclearDiatom" || species_type == "HeteronuclearDiatom" ||
                                  species_type == "LinearTriatomInversionCenter"))
                            {
                                Log<Message>::Error(
                                    "Invalid J entry ", J_str,
                                    " in compound rotational state. LoKI-B only supports compound "
                                    "rotational states for species types with a single rotational quanta.");
                            }

                            J_vals.insert(std::stoi(J_str));
                        }
                        if (J_vals.size() != rot_cnf.size())
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
                if (vib_cnf.size() == 0)
                {
                    throw std::runtime_error("At least one vibrational state is expected by LoKI-B.");
                }

                // Special case treatment for species types with a single
                // vibrational quanta. Compound vibrational states can then be
                // written as a range (e.g. "0-5").
                if (species_type == "HomonuclearDiatom" || species_type == "HeteronuclearDiatom") {
                    std::set<unsigned> v_vals;
                    for (const auto &ventry : vib_cnf)
                    {
                        if (ventry.contains("rotational"))
                        {
                            throw std::runtime_error("Rotational states identifiers are not allowed when "
                                                     "multiple vibrational states are specified.");
                        }
                        const auto &v_str = ventry.at("summary").get<std::string>();
                        v_vals.insert(std::stoi(v_str));
                    }
                
                    // For "v" and "J" we assume that the entries form a continuous value-range.
                    if (v_vals.size() != vib_cnf.size())
                    {
                        throw std::runtime_error("Duplicate v entries encountered.");
                    }
                    auto nv = *v_vals.rbegin() + 1 - *v_vals.begin();
                    if (nv != v_vals.size())
                    {
                        throw std::runtime_error("Expected a contiguous v-range.");
                    }
                    v = std::to_string(*v_vals.begin()) + '-' + std::to_string(*v_vals.rbegin());
                } else {
                    for (nlohmann::json::size_type i = 0; i < vib_cnf.size(); i++) {
                        v.append(vib_cnf[i].at("summary").get<std::string>());
                        if (i < vib_cnf.size() - 1) v.append("|");
                    }
                }
            }
        }
    }
    StateType stateType = J.empty() == false   ? rotational
                          : v.empty() == false ? vibrational
                          : e.empty() == false ? electronic
                                               : charge;
    return StateEntry{id, stateType, gas_name, charge_str, e, v, J};
}

StateEntry propertyStateFromString(const std::string &propertyString)
{
    static const std::regex reState(
        R"(([A-Za-z][A-Za-z0-9]*)\(([-\+]*)(?:\s*,)?\s*([-\+'\[\]/\w\*_|\^]*)\s*(?:,\s*v\s*=\s*([-\+\w\*]+))?\s*(?:,\s*J\s*=\s*([-\+\d\*]+))?\s*\))");
    std::smatch m;

    if (!std::regex_match(propertyString, m, reState))
        return {};

    if (m.str(1).empty())
        return {};

    StateType stateType;

    if (m.str(3).empty())
    {
        stateType = charge;
    }
    else if (m.str(4).empty())
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
