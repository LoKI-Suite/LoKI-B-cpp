/** \file
 *
 *  Unit tests for Range.h
 *
 *  \author Jan van Dijk
 *  \date   May 2024
 */

#include "LoKI-B/Range.h"
#include "tests/TestUtilities.h"
#include <cmath>
#include <memory>

void test1()
{
    using Range = loki::Range;
    using json_type = loki::json_type;

    // value
    const json_type valrange_cnf{ { "value", 42 } };
    std::unique_ptr<const Range> valrange{Range::create(valrange_cnf)};
    test_expr(valrange->size()==1 && valrange->value(0)==42.);

    // linspan
    const json_type linspace_cnf{ {"range","linspace" },{ "min", -1},{ "max", 1 },{ "nvalues", 3 } };
    std::unique_ptr<const Range> linspace{Range::create(linspace_cnf)};
    test_expr(linspace->size()==3);
    for (Range::size_type i=0; i!=linspace->size(); ++i)
    {
        test_expr( std::abs(linspace->value(i)-(int(i)-1)) < 10*std::numeric_limits<double>::epsilon());
    }

    // logspace
    const json_type logspace_cnf{ {"range","logspace" },{ "min", -1},{ "max", 1 },{ "nvalues", 3 } };
    std::unique_ptr<const Range> logspace{Range::create(logspace_cnf)};
    test_expr(logspace->size()==3);
    for (Range::size_type i=0; i!=logspace->size(); ++i)
    {
        test_expr( std::abs(logspace->value(i)-std::pow(10,int(i)-1)) < 10*std::numeric_limits<double>::epsilon());
    }

    // user-specified array
    const json_type array_range_cnf{ -1, 0, 1 };
    std::unique_ptr<const Range> array_range{Range::create(array_range_cnf)};
    test_expr(array_range->size()==3);
    for (Range::size_type i=0; i!=array_range->size(); ++i)
    {
        test_expr( std::abs(array_range->value(i)-(int(i)-1)) < 10*std::numeric_limits<double>::epsilon());
    }
}

int main()
{
    test1();

    test_report;
    return nerrors;
}
