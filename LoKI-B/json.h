#ifndef LOKI_CPP_JSON_H
#define LOKI_CPP_JSON_H

#include "nlohmann/json.hpp"
#include "LoKI-B/Exports.h"
#include <iostream>
#include <string>

namespace loki
{

using json_type = nlohmann::json;

lokib_export json_type read_json_from_stream(std::istream &is);
lokib_export json_type read_json_from_file(const std::string &fname);

} // namespace loki

#endif // LOKI_CPP_JSON_H
