/** \file
 *
 *  Unit tests for invalid inputs.
 *
 *  \author Daan Boer
 *  \date   June 2024
 */

#include "LoKI-B/Environment.h"
#include "LoKI-B/LegacyToJSON.h"
#include "LoKI-B/Simulation.h"
#include "tests/TestUtilities.h"
#include <fstream>
#include <iostream>
#include <stdexcept>

void disable_output()
{
    std::cout.setstate(std::ios_base::failbit);
    std::cerr.setstate(std::ios_base::failbit);
}

void restore_output()
{
    std::cout.clear();
    std::cerr.clear();
}

void should_throw(const std::filesystem::path &test_path, const std::string &expected_error)
{
    std::cout << test_path << std::endl;
    disable_output();

    const auto input = loki::legacyToJSON(test_path);
    try
    {
        loki::Simulation test_sim(test_path, input);
        test_sim.run();
    }
    catch (const std::exception &error)
    {
        restore_output();
        test_expr(error.what() == expected_error);
        return;
    }

    restore_output();
    test_expr(false);
}

void test1()
{
    const std::string fname = loki::getEnvironmentVariable("LOKI_TEST_INPUT_DIR") + "/invalid-inputs/incomplete-populations/input.in";
    std::cout << fname << std::endl;
    std::filesystem::path input_path("/home/daan/git/loki-b/tests/invalid-inputs/incomplete-populations/input.in");
    std::string error("N2(X)\n");
    should_throw(input_path, error);
}

int main()
{
    test1();

    test_report;
    return nerrors;
}
