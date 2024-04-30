#ifndef LOKI_CPP_JSON_H
#define LOKI_CPP_JSON_H

/** Instruct nlohmann::json objects to keep track of their parents, so
 *  sensible error messages can be produced when (for example) access
 *  to an element of a type conversion fails.
 *
 *  See https://json.nlohmann.me/api/macros/json_diagnostics/
 *
 *  gcc is known to produce false positives when the -Warray-bounds
 *  flag is used --- and indeed it does. That flag has been disabled
 *  in CMakeLists.txt for this reason. For more details, see
 *  https://github.com/nlohmann/json/issues/3808
 */
#define JSON_DIAGNOSTICS 1

#include <nlohmann/json.hpp>
#include "LoKI-B/Exports.h"
#include <filesystem>
#include <iostream>
#include <string>

namespace loki
{

using json_type = nlohmann::json;

lokib_export json_type read_json_from_stream(std::istream &is);
lokib_export json_type read_json_from_file(const std::filesystem::path &fname);

} // namespace loki

#endif // LOKI_CPP_JSON_H
