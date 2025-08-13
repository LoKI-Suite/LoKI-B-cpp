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
        return error.what() == expected_error;
    }

    return false;
}

int main(int argc, char **argv)
{
    // Tests that should pass.
    test_expr(test_valid("N2()", R"json({
            "detailed": {
              "type": "HomonuclearDiatom",
              "charge": 0
            },
            "serialized": {
              "composition": { "summary": "N2" }
            }
          })json"_json));
    test_expr(test_valid("N2(+)", R"json({
            "detailed": {
              "type": "HomonuclearDiatom",
              "charge": 1
            },
            "serialized": {
              "composition": { "summary": "N2" }
            }
          })json"_json));
    test_expr(test_valid("N2(++)", R"json({
            "detailed": {
              "type": "HomonuclearDiatom",
              "charge": 2
            },
            "serialized": {
              "composition": { "summary": "N2" }
            }
          })json"_json));
    test_expr(test_valid("N2(-)", R"json({
            "detailed": {
              "type": "HomonuclearDiatom",
              "charge": -1
            },
            "serialized": {
              "composition": { "summary": "N2" }
            }
          })json"_json));
    test_expr(test_valid("N2(--)", R"json({
            "detailed": {
              "type": "HomonuclearDiatom",
              "charge": -2
            },
            "serialized": {
              "composition": { "summary": "N2" }
            }
          })json"_json));
    test_expr(test_valid("N2(X)", R"json({
            "detailed": {
              "type": "HomonuclearDiatom",
              "charge": 0
            },
            "serialized": {
              "composition": { "summary": "N2" },
              "electronic": { "summary": "X" }
            }
          })json"_json));
    test_expr(test_valid("Ar(^1S_0)", R"json({
            "detailed": {
              "type": "AtomLS",
              "charge": 0
            },
            "serialized": {
              "composition": { "summary": "Ar" },
              "electronic": { "summary": "^1S_0" }
            }
          })json"_json));
    test_expr(test_valid("N2(X|A)", R"json({
           "detailed": {
             "type": "HomonuclearDiatom",
             "charge": 0
           },
           "serialized": {
             "composition": { "summary": "N2" },
             "electronic": [
               { "summary": "X" },
               { "summary": "A" }
             ]
           }
         })json"_json));
    test_expr(test_valid("N2(X,v=0)", R"json({
            "detailed": {
              "type": "HomonuclearDiatom",
              "charge": 0
            },
            "serialized": {
              "composition": { "summary": "N2" },
              "electronic": {
                "summary": "X",
                "vibrational": { "summary": "0" }
              }
            }
          })json"_json));
    test_expr(test_valid("N2(X,v=0-4)", R"json({
            "detailed": {
              "type": "HomonuclearDiatom",
              "charge": 0
            },
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
            "detailed": {
              "type": "HomonuclearDiatom",
              "charge": 0
            },
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
            "detailed": {
              "type": "HomonuclearDiatom",
              "charge": 0
            },
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
                             "detailed": {
                               "type": "LinearTriatomInversionCenter",
                               "charge": 0
                             },
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
    test_expr(test_valid("CO2(X,v=010|011)",
                         R"json({
                             "detailed": {
                               "type": "LinearTriatomInversionCenter",
                               "charge": 0
                             },
                             "serialized": {
                               "composition": { "summary": "CO2" },
                               "electronic": {
                                 "summary": "X",
                                 "vibrational": [
                                   { "summary": "010" },
                                   { "summary": "011" }
                                 ]
                               }
                             }
                           })json"_json));

    // Tests that should fail.
    test_expr(
        test_invalid("Rotational states identifiers are not allowed when multiple vibrational states are specified.",
                     R"json({
                       "detailed": {
                         "type": "HomonuclearDiatom",
                         "charge": 0
                       },
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
                             "detailed": {
                               "type": "HomonuclearDiatom",
                               "charge": 0
                             },
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
                             "detailed": {
                               "type": "HomonuclearDiatom",
                               "charge": 0
                             },
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
                             "detailed": {
                               "type": "HomonuclearDiatom",
                               "charge": 0
                             },
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
    test_expr(test_invalid("Expected a contiguous J-range.",
                           R"json({
                             "detailed": {
                               "type": "HomonuclearDiatom",
                               "charge": 0
                             },
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
                             "detailed": {
                               "type": "HomonuclearDiatom",
                               "charge": 0
                             },
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
                             "detailed": {
                               "type": "HomonuclearDiatom",
                               "charge": 0
                             },
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
                             "detailed": {
                               "type": "HomonuclearDiatom",
                               "charge": 0
                             },
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
    test_expr(test_invalid("Invalid J entry 010 in compound rotational state. LoKI-B only supports compound "
                           "rotational states for species types with a single rotational quanta.\n",
                           R"json({
                             "detailed": {
                               "type": "TriatomC2v",
                               "charge": 0
                             },
                             "serialized": {
                               "composition": { "summary": "H2O" },
                               "electronic": {
                                 "summary": "X",
                                 "vibrational": {
                                   "summary": "010",
                                   "rotational": [{ "summary": "010" }, { "summary": "100" }]
                                 }
                                 
                               }
                             }
                           })json"_json));

    test_report;
    return nerrors;
}
