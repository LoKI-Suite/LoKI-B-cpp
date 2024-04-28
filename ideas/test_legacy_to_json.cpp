/** \file
 *
 *  Test the lokib::Event<> template.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2022 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
 *  M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
 *  L. L. Alves
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  \author Jan van Dijk
 *  \date   11 September 2022
 */

#include "ideas/LegacyToJSON.h"
#include "tests/TestUtilities.h"
#include <iostream>
#include <sstream>

void do_test_error(const std::string& str)
{
    ++ntests;
    std::stringstream ss;
    ss << str << '\n';
    try {
        const loki::json_type json = loki::legacyToJSON(ss);
        // if we did *not* get an exception, the test fails.
        ++nerrors;
        std::cout << "ERROR: code was expected to fail. Parsed:\n" << json.dump(2) << std::endl;
    }
    catch(std::exception& exc)
    {
        std::cout << "OK, expected exception: " << exc.what() << std::endl;
    }
}

void test_errors()
{
    do_test_error(":");                   // empty key
    do_test_error(": B");                 // empty key
    do_test_error("  A: B\nC: D");        // bad indentation
    do_test_error("A: B\n  C: D\n E: F"); // bad indentation
    do_test_error("\tA: B");              // tab character in leading whitespace
    do_test_error("A B");                 // expect key: value
    do_test_error("A: B\nA: C");          // duplicate label
    do_test_error("A: B\nA: C");          // duplicate label
    do_test_error("A:  \nA: C");          // duplicate label
    do_test_error("A: B\nA:");            // duplicate label
    do_test_error("A:\n - B\n C: D");     // mixed array/object
    do_test_error("A:\n B: C\n - D");     // mixed array/object
}

void test_correctness_1()
{
    const std::string str = R"(
root:
 S: D
 Bf: false
 Bt: true
 I: 42
 D: 42.0
 section1:
   - 1
   - 2
   - 3
 section2:
   A: B
   C: D
	)";
    std::stringstream ss;
    ss << str << '\n';
    try {
        const loki::json_type json = loki::legacyToJSON(ss);
        std::cout << "OK. Parsed:\n" << json.dump(2) << std::endl;
	test_expr(json["root"]["S"]=="D");
	test_expr(json["root"]["Bf"]==false);
	test_expr(json["root"]["Bt"]==true);
	test_expr(json["root"]["I"]==42);
	test_expr(json["root"]["D"]==42.0);
	test_expr(json["root"]["section1"].size()==3);
	test_expr(json["root"]["section1"][0]==1);
	test_expr(json["root"]["section1"][1]==2);
	test_expr(json["root"]["section1"][2]==3);
	test_expr(json["root"]["section2"].size()==2);
	test_expr(json["root"]["section2"]["A"]=="B");
	test_expr(json["root"]["section2"]["C"]=="D");
    }
    catch(std::exception& exc)
    {
        std::cout << "ERROR: unexpected exception: " << exc.what() << std::endl;
        ++nerrors;
    }
}

int main()
{
    test_errors();
    test_correctness_1();

    test_report;
    return nerrors;
}
