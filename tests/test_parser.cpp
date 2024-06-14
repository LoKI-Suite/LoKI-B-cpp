/** \file
 *  Unit tests for Parse.h
 *  \todo Add more tests, use a testing framework that can be reused for other
 *  test file as well
 *
 *  \author Jan van Dijk
 *  \date   November 2020
 */

#include "LoKI-B/Parse.h"
#include "tests/TestUtilities.h"
#include <iostream>
#include <sstream>

inline std::string escapeNewlineTab(const std::string& src)
{
    std::string copy(src);
    loki::Parse::searchAndReplaceInPlace(copy,"\n","\\n");
    loki::Parse::searchAndReplaceInPlace(copy,"\t","\\t");
    return copy;
}


void testRemoveComments(const std::string& src, const std::string& expected)
{
    std::stringstream ss{src};
    std::string res;
    loki::Parse::removeComments(ss,res);
    test_expr(res==expected);
    std::cout << " Source string (escaped) is '" << escapeNewlineTab(src)
        << "'. Result (escaped) is '" << escapeNewlineTab(res)
        << "', expected '" << escapeNewlineTab(expected) << "'" << std::endl;
}

void testFunctionCall(const std::string& str, const std::string& name, const std::vector<loki::json_type>& args)
try {
    std::cout << "Testing: " << str << std::endl;
    loki::Parse::FunctionCall res(args.size(),str);
    test_expr(res.name() == name);
    test_expr(res.args() == args);
}
catch(std::exception& exc)
{
    std::cerr << "ERROR: " << exc.what() << std::endl;
    ++nerrors;
}

void testFunctionCallException(const std::string& str, unsigned arity=0)
try {
    loki::Parse::FunctionCall res(arity,str);
    // an exception should have been thrown. If we arrive here, the test failed.
    std::cerr << "ERROR: function call '" << str << "' unexpectedly succeeded." << std::endl;
    ++nerrors;
}
catch(std::exception& exc)
{
    std::cout << "OK (expected failure): '" << str << "': '" << exc.what() << "'." << std::endl;
    ++ntests;
}

void testFunctionCalls()
{
    testFunctionCall("f()", "f", {} );
    testFunctionCall("_f()", "_f", {} );
    testFunctionCall("_1f()", "_1f", {} );
    testFunctionCall("f(1)", "f", { 1 } );
    testFunctionCall("f(-1.2e3)", "f", { -1.2e3 } );
    testFunctionCall("f(-1.2e-3)", "f", { -1.2e-3 } );
    testFunctionCall("f(x)", "f", { "x" } );
    testFunctionCall("f(1,2)", "f", { 1,2 } );
    testFunctionCall("f(x,y)", "f", { "x","y" } );
    testFunctionCall("f(x,1)", "f", { "x",1 } );
    testFunctionCall("f(x,y,z)", "f", { "x","y","z" } );
    testFunctionCall("f(1,2,3)", "f", { 1,2,3 } );
    testFunctionCall("f(-1,-2,-3)", "f", { -1,-2,-3 } );
    testFunctionCall("f(x,_y,3)", "f", { "x","_y",3 } );
    testFunctionCall("f(aap,noot,mies,wim)", "f", { "aap","noot","mies", "wim" } );

    testFunctionCallException("");
    testFunctionCallException(" ");
    testFunctionCallException("f");
    testFunctionCallException("%()");
    testFunctionCallException("1()");
    testFunctionCallException("f(");
    testFunctionCallException("f)");
    testFunctionCallException("f(");
    testFunctionCallException("f(%)",1);
    testFunctionCallException("f(1a)",1);
    testFunctionCallException("f(1 2)",2);
    testFunctionCallException("f(1,g())",2);
    // syntax OK, but wrong arity:
    testFunctionCallException("f()",1);
    testFunctionCallException("f(1)",2);
    testFunctionCallException("f(1,2)",3);
    testFunctionCallException("f(1,2,3)",4);

    test_expr( (
        loki::Parse::makeJsonFromFunctionCall("f()", "function", { })
        ==
        loki::json_type{ {"function","f" } }
    ) );
    test_expr( (
        loki::Parse::makeJsonFromFunctionCall("sin(omegat)", "function", { "phase" })
        ==
        loki::json_type{ {"function","sin" }, {"phase", "omegat"} }
    ) );
    test_expr( (
        loki::Parse::makeJsonFromFunctionCall("atan2(4,3)", "function", { "y", "x" })
        ==
        loki::json_type{ {"function","atan2" }, {"y", 4}, {"x", 3} }
    ) );
    test_expr( (
        loki::Parse::makeJsonFromFunctionCall("logspan(-3,3,7)", "range", { "min", "max", "nvalues" })
        ==
        loki::json_type{ {"range","logspan" },{ "min", -3},{ "max", 3 },{ "nvalues", 7 } }
    ) );
}

int main()
{
    testRemoveComments(" \t A    ", " \t A");
    testRemoveComments(" \t A    \n B     \t ", " \t A\n B");
    testRemoveComments("\n% comment \n\n", "");
    testRemoveComments("\n \t% comment \n\n", "");
    testRemoveComments("A % comment\nB", "A\nB");
    testFunctionCalls();

    test_report;

    return nerrors;
}
