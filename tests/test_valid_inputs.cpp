/** \file
 *
 *  Unit and regression tests for valid inputs.
 *
 *  This test will perform all tests defined in the "valid-inputs" directory. A
 *  correct test is defined by a folder that contains an input file "input.in"
 *  and an expected output file "output.json". The test will run a simulation
 *  based on "input.in", and compare its output to "output.json". For now this
 *  is a simple equivalency check on the resulting JSON output object and the
 *  expected JSON object defined in "output.json".
 *
 *  \author Daan Boer
 *  \date   December 2024
 */

#include "LoKI-B/Environment.h"
#include "LoKI-B/LegacyToJSON.h"
#include "LoKI-B/Output.h"
#include "LoKI-B/Simulation.h"
#include "LoKI-B/json.h"
#include "tests/TestUtilities.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
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

bool compare_floats(double a, double b, double rel_tol, double abs_tol)
{
    const auto equal = std::fabs(a - b) <= std::max(rel_tol * std::max(std::fabs(a), std::fabs(b)), abs_tol);
    if (!equal) {
        std::cerr << "Floats " << a << " and " << b << " are not equal." << std::endl;
    }
    return equal;
}

bool json_equal(const nlohmann::json &j1, const nlohmann::json &j2, double rel_tol = 1e-9, double abs_tol = 1e-12)
{
    if (j1.type() != j2.type())
    {
        return false;
    }

    switch (j1.type())
    {
    case nlohmann::json::value_t::object:
        if (j1.size() != j2.size())
        {
            return false;
        }
        for (const auto &item : j1.items())
        {
            std::cerr << "Traversing key " << item.key() << std::endl;
            if (!j2.contains(item.key()) || !json_equal(item.value(), j2.at(item.key()), rel_tol, abs_tol))
            {
                return false;
            }
        }
        return true;

    case nlohmann::json::value_t::array:
        if (j1.size() != j2.size())
        {
            return false;
        }
        for (size_t i = 0; i < j1.size(); ++i)
        {
            std::cerr << "At index " << i << std::endl;
            if (!json_equal(j1[i], j2[i], rel_tol, abs_tol))
            {
                return false;
            }
        }
        return true;

    case nlohmann::json::value_t::number_float:
        return compare_floats(j1.get<double>(), j2.get<double>(), rel_tol, abs_tol);

    default:
        return j1 == j2;
    }
}

void execute_test(const fs::path &test_dir)
{
    const auto input_path = test_dir / "input.in";
    const auto expected_output = loki::read_json_from_file(test_dir / "output.json");

    const auto test_output = std::make_unique<loki::json_type>();

    disable_output();

    const auto input = loki::legacyToJSON(input_path);
    try
    {
        loki::Simulation test_sim(input_path, input);
        test_sim.obtainedResults().addListener(
            &loki::Output::saveCycle, new loki::JsonOutput(*test_output, input, &test_sim.workingConditions()));
        test_sim.run();
    }
    catch (const std::exception &error)
    {
        restore_output();
        std::cerr << error.what() << std::endl;
        test_expr(false);
        return;
    }

    restore_output();
    test_expr(json_equal(expected_output, *test_output));
}

void regenerate_output(const fs::path &test_dir)
{
    const auto input_path = test_dir / "input.in";
    const auto test_output = std::make_unique<loki::json_type>();

    disable_output();

    const auto input = loki::legacyToJSON(input_path);
    loki::Simulation test_sim(input_path, input);
    test_sim.obtainedResults().addListener(&loki::Output::saveCycle,
                                           new loki::JsonOutput(*test_output, input, &test_sim.workingConditions()));
    test_sim.run();

    restore_output();

    std::ofstream ofs(test_dir / "output.json");
    ofs << test_output->dump(2);
}

int main()
{
    const fs::path base_dir = loki::getEnvironmentVariable("LOKI_TEST_INPUT_DIR");

    unsigned counter = 1;

    for (const auto &dir : fs::directory_iterator(base_dir / "valid-inputs"))
    {
        if (dir.is_directory() && fs::exists(dir.path() / "input.in") && fs::exists(dir.path() / "output.json"))
        {
            std::cout << counter++ << ". " << dir.path().filename().c_str() << std::endl;
            execute_test(dir.path());
        }
    }

    test_report;
    return nerrors;
}
