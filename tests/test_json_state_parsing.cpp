#include "LoKI-B/StateEntry.h"
#include "tests/TestUtilities.h"
#include <exception>
#include <sstream>
#include <string_view>

bool test_valid(const std::string_view state_str, const nlohmann::json &state)
{
    const auto entry = loki::entryFromJSON("TestID", state);

    std::ostringstream ss;
    ss << entry;

    return ss.str() == state_str;
}

bool test_invalid(const std::string_view expected_error, const nlohmann::json &state)
{
    try
    {
        const auto entry = loki::entryFromJSON("TestID", state);
        std::cout << entry << std::endl;
    }
    catch (const std::exception &error)
    {
        std::cout << error.what() << std::endl;
        return error.what() == expected_error;
    }

    return false;
}

int main(int argc, char **argv)
{
    // Tests that should pass.
    test_expr(test_valid("N2(X)", R"json({
            "serialized": {
              "composition": { "summary": "N2" },
              "electronic": { "summary": "X" }
            }
          })json"_json));
    test_expr(test_valid("Ar(^1S_0)", R"json({
            "serialized": {
              "composition": { "summary": "Ar" },
              "electronic": { "summary": "^1S_0" }
            }
          })json"_json));
    test_expr(test_valid("N2(X,v=0)", R"json({
            "serialized": {
              "composition": { "summary": "N2" },
              "electronic": {
                "summary": "X",
                "vibrational": { "summary": "0" }
              }
            }
          })json"_json));
    test_expr(test_valid("N2(X,v=0-4)", R"json({
            "serialized": {
              "composition": { "summary": "N2" },
              "electronic": {
                "summary": "X",
                "vibrational": [
                  { "summary": "0" },
                  { "summary": "1" },
                  { "summary": "2" },
                  { "summary": "3" },
                  { "summary": "4" }
                ]
              }
            }
          })json"_json));
    test_expr(test_valid("N2(X,v=0,J=0)", R"json({
            "serialized": {
              "composition": { "summary": "N2" },
              "electronic": {
                "summary": "X",
                "vibrational": {
                  "summary": "0",
                  "rotational": {
                    "summary": "0"
                  }
                }
              }
            }
          })json"_json));
    test_expr(test_valid("N2(X,v=0,J=3-6)", R"json({
            "serialized": {
              "composition": { "summary": "N2" },
              "electronic": {
                "summary": "X",
                "vibrational": {
                  "summary": "0",
                  "rotational": [
                    { "summary": "3" },
                    { "summary": "4" },
                    { "summary": "5" },
                    { "summary": "6" }
                  ]
                }
              }
            }
          })json"_json));
    test_expr(test_valid("CO2(X,v=010,J=0)",
                         R"json({
                             "serialized": {
                               "composition": { "summary": "CO2" },
                               "electronic": {
                                 "summary": "X",
                                 "vibrational": {
                                   "summary": "010",
                                   "rotational": { "summary": "0" }
                                 }
                               }
                             }
                           })json"_json));

    // Tests that should fail.
    test_expr(test_invalid("Exactly one electronic state is expected by LoKI-B.",
                           R"json({
                             "serialized": {
                               "composition": { "summary": "N2" },
                               "electronic": [
                                 { "summary": "X" },
                                 { "summary": "A" }
                               ]
                             }
                           })json"_json));
    test_expr(
        test_invalid("Rotational states identifiers are not allowed when multiple vibrational states are specified.",
                     R"json({
                       "serialized": {
                         "composition": { "summary": "N2" },
                         "electronic": {
                           "summary": "X",
                           "vibrational": [
                             { "summary": "0", "rotational": { "summary": "0" } },
                             { "summary": "1" }
                           ]
                         }
                       }
                     })json"_json));
    test_expr(test_invalid("At least one vibrational state is expected by LoKI-B.",
                           R"json({
                             "serialized": {
                               "composition": { "summary": "N2" },
                               "electronic": {
                                 "summary": "X",
                                 "vibrational": []
                               }
                             }
                           })json"_json));
    test_expr(test_invalid("Expected a contiguous v-range.",
                           R"json({
                             "serialized": {
                               "composition": { "summary": "N2" },
                               "electronic": {
                                 "summary": "X",
                                 "vibrational": [
                                   { "summary": "0" },
                                   { "summary": "2" }
                                 ]
                               }
                             }
                           })json"_json));
    test_expr(test_invalid("Duplicate v entries encountered.",
                           R"json({
                             "serialized": {
                               "composition": { "summary": "N2" },
                               "electronic": {
                                 "summary": "X",
                                 "vibrational": [
                                   { "summary": "0" },
                                   { "summary": "0" }
                                 ]
                               }
                             }
                           })json"_json));
    // FIXME: This test currently does not behave as it should. Parsing of this
    // definition passes and produces "CO2(X,v=10-11)". It should instead throw
    // an error. Should we check for leading zeros?
    // test_expr(test_invalid("",
    //                        R"json({
    //                          "serialized": {
    //                            "composition": { "summary": "CO2" },
    //                            "electronic": {
    //                              "summary": "X",
    //                              "vibrational": [
    //                                { "summary": "010" },
    //                                { "summary": "011" }
    //                              ]
    //                            }
    //                          }
    //                        })json"_json));
    test_expr(test_invalid("Expected a contiguous J-range.",
                           R"json({
                             "serialized": {
                               "composition": { "summary": "N2" },
                               "electronic": {
                                 "summary": "X",
                                 "vibrational": {
                                   "summary": "0",
                                   "rotational": [
                                     { "summary": "0" },
                                     { "summary": "2" }
                                   ]
                                 }
                               }
                             }
                           })json"_json));
    test_expr(test_invalid("Expected a contiguous J-range.",
                           R"json({
                             "serialized": {
                               "composition": { "summary": "N2" },
                               "electronic": {
                                 "summary": "X",
                                 "vibrational": {
                                   "summary": "0",
                                   "rotational": [
                                     { "summary": "0" },
                                     { "summary": "2" }
                                   ]
                                 }
                               }
                             }
                           })json"_json));
    test_expr(test_invalid("Duplicate J entries encountered.",
                           R"json({
                             "serialized": {
                               "composition": { "summary": "N2" },
                               "electronic": {
                                 "summary": "X",
                                 "vibrational": {
                                   "summary": "0",
                                   "rotational": [
                                     { "summary": "0" },
                                     { "summary": "0" }
                                   ]
                                 }
                               }
                             }
                           })json"_json));
    test_expr(test_invalid("At least one rotational state is expected by LoKI-B.",
                           R"json({
                             "serialized": {
                               "composition": { "summary": "N2" },
                               "electronic": {
                                 "summary": "X",
                                 "vibrational": {
                                   "summary": "0",
                                   "rotational": []
                                 }
                               }
                             }
                           })json"_json));

    test_report;
    return nerrors;
}
