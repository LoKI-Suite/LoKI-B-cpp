#include "LoKI-B/Parse.h"

namespace {
    std::string group(const std::string& term) { return "(" + term + ")"; }
}

/** \todo Populate the enries and StoiCoeff vectors with the parse results
 */
void entriesFromStringNew(const std::string stateString, std::vector<loki::StateEntry> entries, std::vector<uint16_t>* stoiCoeff)
{
    /// \todo Remove when integrated:
    using namespace loki;

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
            Enumeration::StateType stateType
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

bool test_parser(const std::string stateString, bool new_parser)
{
    try {
        using namespace loki;
        std::vector<StateEntry> entries;
        std::vector<uint16_t> stoiCoeff;
        std::cout << "Testing: '" << stateString << "'." << std::endl;
        if (new_parser)
        {
            entriesFromStringNew(stateString, entries, &stoiCoeff);
        }
        else
        {
            Parse::entriesFromString(stateString, entries, &stoiCoeff);
        }
        assert(entries.size()==stoiCoeff.size());
        for (std::size_t ndx=0; ndx!=entries.size(); ++ndx)
        {
            std::cout << stoiCoeff[ndx] << '\t' << entries[ndx] << ", type = " << entries[ndx].level << std::endl;
        }
        return true;
    }
    catch (std::exception& exc)
    {
        std::cout << "Error: " << exc.what() << std::endl;
        return false;
    }
}

unsigned nerrors_old=0;
unsigned nerrors_new=0;

void test(const std::string& stateString, bool expected=true)
{
    if (test_parser(stateString,false)!=expected)
    {
        std::cout << "ERROR (OLD PARSER): Unexpected result " << !expected << " for input " << stateString << std::endl;
        ++nerrors_old;
    }
    if (test_parser(stateString,true)!=expected)
    {
        std::cout << "ERROR (NEW PARSER): Unexpected result " << !expected << " for input " << stateString << std::endl;
        ++nerrors_new;
    }
}

int main(int argc, const char* argv[])
{
    // acceptance tests:

    test("N2(++,X,v=0,J=1)",true);
    test("N2(++,X,v=0)",true);
    test("N2(++,X)",true);
    test("N2(X,v=0,J=1)",true);
    test("N2(X,v=0)",true);
    test("N2(X)",true);
    test("e",true);
    test("2 e",true); // NOTE: whitespace between coef and particle is optional
    test("2 N2(X)",true);
    test("2e",true);
    test("2N2(X)",true);

    test("e + N2(X)",true);
    test("N2(X) + N2(+,X)",true);
    test("2 N2(X) + 3 N2(+,X)",true);
    test("2N2(X) + 3N2(+,X)",true);

    // rejection tests:

    // no entries
    test("",false);
    // bad gas name (must start with A-z and consist only of A-z and 0-9)
    test("N_2(X)",false);
    // the same (the charge specifier is part of the state id list at present)
    test("N2+(X)",false);

    // missing or spurious parentheses
    test("N_2(X",false);
    test("N_2(X(",false);
    test("N_2(X())",false);
    test("N_2(X(bla))",false);
    test("N_2(X(bla),v=1)",false);
    // bad charge id:
    test("N2(+-,X,v=0,J=1)",false);
    test("N2(bla,X,v=0,J=1)",false);
    // bad elec id:
    test("N2(++,e=X,v=0,J=1)",false);
    test("N2(--,e=X,v=0,J=1)",false);
    // bad vib id:
    test("N2(++,X,0,J=1)",false);
    test("N2(--,X,0,J=1)",false);
    // bad rot id:
    test("N2(++,X,v=0,1)",false);
    test("N2(--,X,v=0,1)",false);
    // missing elec id:
    test("N2(++,v=0)",false);
    test("N2(--,v=0)",false);
    // missing elec id:
    test("N2(++,v=0,J=1)",false);
    test("N2(--,v=0,J=1)",false);
    // mssing vib id:
    test("N2(++,X,J=1)",false);
    test("N2(--,X,J=1)",false);
    // mssing vib id:
    test("N2(++,X,J=1)",false);
    test("N2(--,X,J=1)",false);
    // wrong order
    test("N2(++,X,J=1,v=0)",false);
    test("N2(--,J=1,X,v=0)",false);
    test("N2(++,J=1,X,v=0)",false);
    test("N2(--,J=1,X,v=0)",false);
    test("N2(++,J=1,v=0,X)",false);
    test("N2(--,J=1,v=0,X)",false);
    test("N2(X,J=1,v=0)",false);
    test("N2(J=1,X,v=0)",false);
    test("N2(J=1,X,v=0)",false);
    test("N2(J=1,X,v=0)",false);
    test("N2(J=1,v=0,X)",false);
    test("N2(J=1,v=0,X)",false);
    // trailing characters
    test("N2(++,X,v=0,J=1,PSV=1)",false);
    test("N2(--,X,v=0,J=1,PSV=1)",false);

    // spurious coef or spurious plus, missing particle:
    test("2 2 N2(X)",false);
    test("N2(X) 2",false);
    test("N2(X) + + N2(X)",false);
    test("+ N2(X)",false);
    test("N2(X)+",false);
    test("N2(X) +",false);

    std::cout << "Number of errors (old parser): " << nerrors_old << std::endl;
    std::cout << "Number of errors (new parser): " << nerrors_new << std::endl;
    return 0;
}
