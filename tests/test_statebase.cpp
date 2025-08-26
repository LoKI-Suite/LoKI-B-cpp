/** \file
 *  Unit tests for class Gas::State.
 *  \todo Add more tests, use a testing framework that can be reused for other
 *  test file as well
 *
 *  \author Jan van Dijk
 *  \date   July 2020
 */

#include "LoKI-B/Gas.h"
#include "LoKI-B/GasProperties.h"
#include "LoKI-B/Log.h"
#include <sstream>

unsigned ntests=0;
unsigned nerrors=0;

void test_state_string(const loki::GasProperties &gasProps, const std::string str, bool should_pass, std::string expected = "")
{
    ++ntests;
    try
    {
        using namespace loki;

        // std::cout << "State string: '" << str << "'." << std::endl;
        const StateEntry e = propertyStateFromString(str);
        // std::cout << "Entry: '" << e << " of type " << e.m_level << "." << std::endl;
        loki::Gas gas(gasProps, e.m_gasName);
        // std::cout << "Gas name: '" << gas.name << "'." << std::endl;
        const Gas::State root{&gas};
        const Gas::State* s{gas.ensureState(e)};

        std::ostringstream entry_stream;
        entry_stream << e;

        std::ostringstream state_stream;
        state_stream << *s;

        if (entry_stream.str() != state_stream.str())
            Log<Message>::Error("Serialized entry and state are not equal.");

        if (state_stream.str() != expected)
            Log<Message>::Error("Serialized entry is not equal to expected result.");

        // std::cout << "State: " << *s << std::endl;
    }
    catch (std::exception &exc)
    {
        if (should_pass)
        {
            std::cerr << " *** ERROR: Unexpected failure for '" << str << "', got: " << exc.what() << std::endl;
            ++nerrors;
        }
        else
        {
            std::cout << " *** OK, passed test for '" << str << "' (expected failure)." << std::endl;
        }
        return;
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
    // We do not care about the gas properties, but mass is mandatory
    loki::GasProperties gasProps;
    gasProps.set("mass", "N2", 4.651834066656000e-26);

    test_state_string(gasProps, "N2", true);   // root N2 state
    test_state_string(gasProps, "N2()", true); // neutral charge state
    test_state_string(gasProps, "N2(X)", true);
    test_state_string(gasProps, "N2(X,v=0)", true);
    test_state_string(gasProps, "N2(X,v=0,J=0)", true);
    test_state_string(gasProps, "N2(+)", true);
    test_state_string(gasProps, "N2(--)", true);

    // TODO: This currently passes, do we want it to fail?
    test_state_string(gasProps, "N2(v=0)", true); // produces state type electronic (not vibrational)

    // Things that should fail:
    test_state_string(gasProps, "N2(", false);       // string is incomplete
    test_state_string(gasProps, "N2)", false);       // string is incomplete
    test_state_string(gasProps, "2N(X)", false);     // gas name should start with a capital
    test_state_string(gasProps, "N2(X,J=0)", false); // missing vibrational specifier

    std::cout << " *** Number of tests: " << ntests << std::endl;
    std::cout << " *** Number of errors: " << nerrors << std::endl;

    return nerrors;
}
