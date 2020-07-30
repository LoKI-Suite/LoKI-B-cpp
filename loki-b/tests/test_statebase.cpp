/** \file
 *  Unit tests for class StateBase.
 *  \todo Add more tests, use a testing framework that can be reused for other
 *  test file as well
 *
 *  \author Jan van Dijk
 *  \date   July 2020
 */

#include "LoKI-B/GasBase.h"
#include "LoKI-B/Parse.h"

// GasBase has a protected destructor. Derive a dummy class to allow the creation of an object.
class TestGas : public loki::GasBase { public: TestGas(const std::string& name) : loki::GasBase(name) {} };

unsigned ntests=0;
unsigned nerrors=0;

void test_state_string(const std::string str, bool should_pass)
{
    ++ntests;
    try {
        using namespace loki;

//        std::cout << "State string: '" << str << "'." << std::endl;
        const StateEntry e = Parse::propertyStateFromString(str);
//        std::cout << "Entry: '" << e << "." << std::endl;
        TestGas gas(e.gasName);
  //      std::cout << "Gas name: '" << gas.name << "'." << std::endl;
        const GasBase::StateBase s{e,e.level,&gas,nullptr};
    //    std::cout << "State:" << std::endl;
      //  std::cout << s << std::endl;
    }
    catch(std::exception& exc)
    {
        if (should_pass)
        {
            std::cerr << " *** ERROR: Unexpected failure for '" << str << "', got: " << exc.what() << std::endl;
            ++nerrors;
        }
        else
        {
            std::cout << " *** OK, passed test for '" << str << "' (expected failuer)."
                << std::endl << "Received: " << exc.what() << std::endl;
        }
    }
    // if we get here, no exception was thrown
    if (should_pass)
    {
        std::cout << " *** OK, passed test for '" << str << "'." << std::endl;
    }
    else
    {
        std::cerr << " *** ERROR: test was expected to fail for '" << str << "'." << std::endl;
        ++nerrors;
    }
}

int main()
{
    test_state_string("N2(X)",true);
    test_state_string("N2(X,v=0)",true);
    test_state_string("N2(X,v=0,J=0)",true);

    // QUESTION: This results in state type 'none'. Should an attempt to create such
    //           state not throw an exception?
    test_state_string("N2",false);

    // Things that should fail:

    test_state_string("N2(",false);   // string is incomplete
    test_state_string("N2)",false);   // string is incomplete
    test_state_string("N2()",false);   // string is incomplete
    test_state_string("2N(X)",false); // 2 is ignored
    test_state_string("N2(v=0)",false); // produces stete type electronic (not vibrational)
    test_state_string("N2(X,J=0)",false); // produces stete type electronic (not rotational).

    std::cout << " *** Number of tests: " << ntests << std::endl;
    std::cout << " *** Number of errors: " << nerrors << std::endl;

    return nerrors;
};
