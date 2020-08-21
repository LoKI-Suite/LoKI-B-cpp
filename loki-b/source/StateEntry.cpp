#include "LoKI-B/StateEntry.h"
#include "LoKI-B/Parse.h"
#include "LoKI-B/StandardPaths.h"
#include <regex>
#include <fstream>

namespace loki {

StateEntry::StateEntry()
    : level(none)
{
}

StateEntry::StateEntry(StateType level, const std::string &gasName, const std::string &charge,
                   const std::string &e, const std::string &v, const std::string &J)
    : level(level),
    charge(charge),
    gasName(gasName),
    e(e),
    v(v),
    J(J)
{
}

/// \todo Update: chare may also be a wildcard (?)
bool StateEntry::hasWildCard()
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
    os << entry.gasName << '(';

    if (!entry.charge.empty()) os << entry.charge << ',';

    os << entry.e << ",v=" << entry.v << ",J=" << entry.J << ')';

    return os;
}

void entriesFromStringOld(const std::string &statesString, std::vector<StateEntry> &entries,
                              std::vector<uint16_t> *stoiCoeff = nullptr)
{
    static const std::regex reState(
        R"((\d*)\s*([A-Za-z][A-Za-z0-9]*)\(([-\+]?)\s*,?\s*([-\+'\[\]/\w]+)\s*(?:,\s*v\s*=\s*([-\+\w]+))?\s*(?:,\s*J\s*=\s*([-\+\d]+))?\s*)");

    std::regex_iterator<std::string::const_iterator> rit(statesString.begin(), statesString.end(), reState);
    std::regex_iterator<std::string::const_iterator> rend;

    if (rit == rend)
    {
        throw std::runtime_error("No entries found.");
    }

    while (rit != rend)
    {
    try {
        StateType stateType;

        if (rit->str(2).empty())
        {
            throw std::runtime_error("Empry gas name.");
        }
        if (rit->str(4).empty())
        {
            throw std::runtime_error("Electron state not specified.");
        }

        if (rit->str(5).empty())
        {
            stateType = electronic;
        }
        else if (rit->str(6).empty())
        {
            stateType = vibrational;
        }
        else
        {
            stateType = rotational;
        }

        if (stoiCoeff != nullptr)
        {
            if (rit->str(1).empty())
            {
                stoiCoeff->emplace_back(1);
            }
            else
            {
                std::stringstream ss(rit->str(1));
                uint16_t coeff;

                ss >> coeff;
                stoiCoeff->emplace_back(coeff);
            }
        }

        entries.emplace_back(stateType, rit->str(2), rit->str(3), rit->str(4), rit->str(5), rit->str(6));
        ++rit;
    }
    catch(std::exception& exc)
    {
        throw std::runtime_error("While parsing particle '" + rit->str(0) + "':\n" + std::string{exc.what()});
    }
    }
}

namespace {

    // Parenthesize a string and return the result. This makes a regex group
    // and is used for readability.
    std::string group(const std::string& term) { return "(" + term + ")"; }

}

void entriesFromStringNew(const std::string stateString, std::vector<StateEntry> entries, std::vector<uint16_t>* stoiCoeff)
{
    std::cout << " * Testing '" << stateString << "'." << std::endl;

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
    try {
        // as a result of the special handling of the electron, two types of
        // submatch-sequences are possible, as shown below. Also the submatch
        // indices are indicated.
        //       0        1       2         3   4      5
        // a. <group> [  <coef> <group> [  'e'  <>     <>    ] ]
        // b. <group> [  <coef> <group> [  <>  <gas> <state> ] ]

        // the stoichiometric coefficient (empty corresponds to 1).
        const std::string c = res[1];
        std::string g,s;
        // state id components:
        std::string q, e, v, J;
        if (res.str(3).empty())
        {
            g = res[4];
            s = res[5];

        try {

            // now parse the state id list
            std::smatch state_res;
            std::string stateRemainder = s;
            const std::regex charge_expr{"^"+group("\\+*|\\-*") + ","};
            if (std::regex_search(stateRemainder, state_res, charge_expr))
            {
                q = state_res[1];
                stateRemainder = state_res.suffix().str();
            }
            const std::regex elec_expr{"^"+group("[^,=]+")+",?"};
            if (std::regex_search(stateRemainder,state_res,elec_expr))
            {
                e = state_res[1];
                stateRemainder = state_res.suffix().str();
            }
            else
            {
                throw std::runtime_error("Bad electronic id, starting at '" + stateRemainder + ".");
            }
            const std::regex vib_expr{"^v="+group("[^,=]+")+",?"};
            if (!stateRemainder.empty())
            {
                if (std::regex_search(stateRemainder,state_res,vib_expr))
                {
                    v = state_res[1];
                    stateRemainder = state_res.suffix().str();
                }
                else
                {
                    throw std::runtime_error("Bad vibrational id, starting at '" + stateRemainder + ".");
                }
            }
            const std::regex rot_expr{"^J="+group("[^,=]+")+",?"};
            if (!stateRemainder.empty())
            {
                if (std::regex_search(stateRemainder,state_res,rot_expr))
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
            StateType stateType
                = J.empty()==false ? rotational
                : v.empty()==false ? vibrational
                : e.empty()==false ? electronic
                : charge;
            entries.push_back(StateEntry(stateType,g,q,e,v,J));
            if (stoiCoeff)
            {
                stoiCoeff->push_back(c.empty() ? 1 : std::stoi(c));
            }
        }
        catch (std::exception& exc)
        {
            throw std::runtime_error("While parsing state identifier list '" + s + "':\n" + std::string(exc.what()));
        }
        }
        else
        {
            g = res[3];
            s = std::string{};
            if (g!="e")
            {
                throw std::logic_error("Expected 'e', found '" + g + ".");
            }
            std::cout << "WARNING: ADDING AN ELECTRON TO THE STATE ENTRY LIST." << std::endl;
            entries.push_back(StateEntry(charge,g,"-",std::string{},std::string{},std::string{}));
            if (stoiCoeff)
            {
                stoiCoeff->push_back(c.empty() ? 1 : std::stoi(c));
            }
        }
        std::cout << "stoich = '" << c
            << "', particle = '" << g << "'."
            << "', state id = '" << s
                << "', q = '" << q
                << "', e = '" << e
                << "', v = '" << v
                << "', J = '" << J
                << "'" << std::endl;
        remainder = res.suffix().str();
    }
    catch(std::exception& exc)
    {
        throw std::runtime_error("While parsing particle '" + part_remainder + "':\n" + std::string{exc.what()});
    }
    }
    // incomplete parse?
    if (!remainder.empty())
    {
        // if parsing fails at the beginning, show the original string without the
        // artificial '+ ' prepended.
        throw std::runtime_error("Parsing of stoichiometric array failed at '"
            + (remainder==parseString ? stateString : remainder) + "'.");
    }
}

StateEntry entryFromJSON(const json_type& cnf)
{
        const std::string gasName = cnf.at("particle");
        const int charge_int = cnf.at("charge").get<int>();
        const std::string charge_str = charge_int ? std::to_string(charge_int) : std::string{};
        const json_type& descr = cnf.at("descriptor");
        if (descr.contains("states"))
        {
            // We have an array of state objects, instead of a single one.
            // loki-b expects a single string for "e", "v" and "J" (the latter
            // can be empty. We need to 'pessimize' the input by canoncatenating
            // the info into single string. The question here is what loki-b
            // supports. Is it allowed, for example, that the "e" fields are
            // different? For now we assume that everything is possible and we
            // simply concatenate all unique "e"'s while for "v" and "J" we assume
            // that the entries form a continuous value-range.
            std::set<std::string> e_vals;
            std::set<unsigned> v_vals;
            std::set<unsigned> J_vals;
            for (json_type::const_iterator s = descr.at("states").begin(); s!= descr.at("states").end(); ++s)
            {
                e_vals.insert(s->at("e").get<std::string>());
                if (s->contains("v"))
                {
                    v_vals.insert(s->at("v").get<int>());
                }
                if (s->contains("J"))
                {
                    J_vals.insert(s->at("J").get<int>());
                }
            }
            std::string e;
            if (e_vals.size()!=1)
            {
                throw std::runtime_error("Expected a unique electronic state identifier.");
            }
            else
            {
                e = *e_vals.begin();
            }
            std::string v;
            if (v_vals.size()==1)
            {
                v = *v_vals.begin();
            }
            else if (v_vals.size()>1)
            {
                int nv = *v_vals.rbegin()+1-*v_vals.begin();
                if (nv!=v_vals.size())
                {
                    throw std::runtime_error("Expected a contiguous v-range.");
                }
                v = std::to_string(*v_vals.begin()) + '-' + std::to_string(*v_vals.rbegin());
            }
            std::string J;
            if (J_vals.size()==1)
            {
                J = *J_vals.begin();
            }
            else if (J_vals.size()>1)
            {
                int nJ = *J_vals.rbegin()+1-*J_vals.begin();
                if (nJ!=J_vals.size())
                {
                    throw std::runtime_error("Expected a contiguous J-range.");
                }
                J = std::to_string(*J_vals.begin()) + '-' + std::to_string(*J_vals.rbegin());
            }
            StateType stateType
                = J.empty()==false ? rotational
                : v.empty()==false ? vibrational
                : e.empty()==false ? electronic
                : charge;
            return StateEntry{stateType,gasName,charge_str,e,v,J};
        }
        else
        {
            const std::string e{descr.at("e").get<std::string>()};
            const std::string v{descr.contains("v") ? (
                descr.at("v").type()==json_type::value_t::string
                    ? descr.at("v").get<std::string>()
                    : std::to_string(descr.at("v").get<int>())
            ) : std::string{} };
            const std::string J{descr.contains("J") ? std::to_string(descr.at("J").get<int>()) : std::string{} };
            /** \todo Check the precise semantics of the next line. Is it possible that J
             *        is specified, but not v? Is e always specified?
             */
            StateType stateType
                = descr.contains("J") ? rotational
                : descr.contains("v") ? vibrational
                : descr.contains("e") ? electronic
                : charge;
            return StateEntry{stateType,gasName,charge_str,e,v,J};
        }
}

bool entriesFromJSON(const json_type& cnf, std::vector<StateEntry> &entries,
                              std::vector<uint16_t> *stoiCoeff)
{
    for (json_type::const_iterator i=cnf.begin(); i!=cnf.end(); ++i)
    {
        if (i->at("particle").get<std::string>() == "e")
        {
            continue;
        }
        /** \todo It appeaes that the present JSON simply repeats the particle
         *        object when it appears more than once. Then the stoichiometric
         *        coefficient of each entry will be one, and we hope that 'e + e'
         *        will be handled the same way as '2 e' (JvD).
         */
        entries.push_back(entryFromJSON(*i));
        if (stoiCoeff)
        {
                stoiCoeff->push_back(1);
        }
    }
    return true;
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
    return {stateType, m.str(1), m.str(2), m.str(3), m.str(4), m.str(5)};
}

bool statePropertyFile(const std::string &fileName, std::vector<std::pair<StateEntry, double>> &entries)
{
    const std::string inputPath = INPUT "/";

    std::ifstream in(inputPath + fileName);

    if (!in.is_open())
        return false;

    std::string line;

    while (std::getline(in, line))
    {
        line = Parse::removeComments(line);

        if (line.size() < 3)
            continue;

        std::string stateString, valueString;

        if (!Parse::stateAndValue(line, stateString, valueString))
            return false;

        double value;

        if (!Parse::getValue(valueString, value))
            return false;

        entries.emplace_back(propertyStateFromString(stateString), value);
    }

    return true;
}


} // namespace loki
