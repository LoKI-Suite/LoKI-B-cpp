/** \file
 *  Unit tests for Parse.h
 *  \todo Add more tests, use a testing framework that can be reused for other
 *  test file as well
 *
 *  \author Jan van Dijk
 *  \date   November 2020
 */

#include "LoKI-B/Parse.h"
#include <iostream>
#include <sstream>

unsigned ntests=0;
unsigned nerrors=0;

inline std::string escapeNewlineTab(const std::string& src)
{
    std::string copy(src);
    loki::Parse::searchAndReplaceInPlace(copy,"\n","\\n");
    loki::Parse::searchAndReplaceInPlace(copy,"\t","\\t");
    return copy;
}


void testRemoveComments(const std::string& src, const std::string& expected)
{
    ++ntests;
    std::stringstream ss{src};
    std::string res;
    loki::Parse::removeComments(ss,res);
    if (res==expected)
    {
        std::cout << " OK, test passed.";
    }
    else
    {
        std::cerr << " ERROR: test failed.";
        ++nerrors;
    }
    std::cout << " Source string (escaped) is '" << escapeNewlineTab(src)
        << "'. Result (escaped) is '" << escapeNewlineTab(res)
        << "', expected '" << escapeNewlineTab(expected) << "'" << std::endl;
}

int main()
{
    testRemoveComments(" \t A    ", " \t A");
    testRemoveComments(" \t A    \n B     \t ", " \t A\n B");
    testRemoveComments("\n% comment \n\n", "");
    testRemoveComments("\n \t% comment \n\n", "");
    testRemoveComments("A % comment\nB", "A\nB");

    std::cout << " *** Number of tests: " << ntests << std::endl;
    std::cout << " *** Number of errors: " << nerrors << std::endl;

    return nerrors;
}
