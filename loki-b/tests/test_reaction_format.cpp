#include "LoKI-B/Parse.h"
#include "LoKI-B/GasBase.h"

bool test_parser(const std::string stateString, bool new_parser)
{
    try {
        using namespace loki;
        std::vector<StateEntry> entries;
        std::vector<uint16_t> stoiCoeff;
        std::cout << "Testing: '" << stateString << "'." << std::endl;
        if (new_parser)
        {
            loki::entriesFromStringNew(stateString, entries, &stoiCoeff);
        }
        else
        {
            loki::entriesFromStringOld(stateString, entries, &stoiCoeff);
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
