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
#include <iterator>
#include <Eigen/Dense>

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

std::string load_error(const std::string &error_file) {
    std::ifstream file(error_file);
    return std::string(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());
}

void should_throw(const fs::path &test_dir)
{
    const auto input_path = test_dir / "input.in";
    const auto expected_error = load_error(test_dir / "error.out");

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

    unsigned counter = 1;

    for (const auto &dir : fs::directory_iterator(base_dir / "invalid-inputs"))
    {
        if (dir.is_directory() && fs::exists(dir.path() / "input.in") && fs::exists(dir.path() / "error.out"))
        {
            std::cout << counter++ << ". " << dir.path().filename().c_str() << std::endl;
            should_throw(dir.path());
        }
    }

    test_report;
    return nerrors;
}
