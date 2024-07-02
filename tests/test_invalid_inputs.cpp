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
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace fs = std::filesystem;

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

void should_throw(const std::string &test_name, const std::string &expected_error)
{
    std::cout << test_name << std::endl;

    fs::path input_path(loki::getEnvironmentVariable("LOKI_TEST_INPUT_DIR"));
    input_path = input_path / "invalid-inputs" / test_name / "input.in";

    disable_output();

    const auto input = loki::legacyToJSON(input_path);
    try
    {
        loki::Simulation test_sim(input_path, input);
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

int main()
{
    const fs::path base_dir = loki::getEnvironmentVariable("LOKI_TEST_INPUT_DIR");

    for (const auto &dir : fs::directory_iterator(base_dir / "invalid-inputs"))
    {
        if (dir.is_directory() && fs::exists(dir.path() / "input.in"))
        {
            should_throw(dir.path().filename(), "N2(X)\n");
        }
    }

    test_report;
    return nerrors;
}
