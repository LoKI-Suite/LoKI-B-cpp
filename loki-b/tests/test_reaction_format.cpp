#include "LoKI-B/Parse.h"

void test(const std::string stateString)
{
    using namespace loki;
    std::vector<StateEntry> entries;
    std::vector<uint16_t> stoiCoeff;
    std::cout << "Testing: '" << stateString << "'." << std::endl;
    bool res = Parse::entriesFromString(stateString, entries, &stoiCoeff);
    if (!res)
    {
       std::cout << "Parsing failed, input was: '" << stateString << "'." << std::endl;
       return;
    }
    assert(entries.size()==stoiCoeff.size());
    for (std::size_t ndx=0; ndx!=entries.size(); ++ndx)
    {
        std::cout << stoiCoeff[ndx] << '\t' << entries[ndx] << ", type = " << entries[ndx].level << std::endl;
    }
}
int main()
{
    test("e + N2(X)");
    test("e + 2N2(X)");
    test("e + 2 N2(X)"); // ! 2 is ignored
    test("e + N2(X"); // ! missing ) is not detected
    test("e + N2(X + Ar(X)"); // ! missing ) is not detected
    test("N2(v=1)"); // missing electronic id, v-spec is ignored
    test("N2(J=1)"); // missing electronic and vibrational id's, J-spec is ignored
    test("N2(v=1,J=1)"); // missing electronic id, v-spec and J-spec are ignored
    test("N2(+,X)");
    test("N2(-,X)");
    test("N2(+)"); // + is interpreted as electronic id
    test("N2(-)"); // - is interpreted as electronic id
    test("N2(++)"); // second + is interpreted as electronic id, missing , is ignored
    test("N2(--)"); // second - is interpreted as electronic id, missing , is ignored
    return 0;
}
