#ifndef LOKI_CPP_JSON_H
#define LOKI_CPP_JSON_H

#include "nlohmann/json.hpp"
#include <iostream>
#include <string>

namespace loki {

using json_type = nlohmann::json;

json_type read_json_from_stream(std::istream& is);
json_type read_json_from_file(const std::string& fname);

}; // namespace loki

#endif // LOKI_CPP_JSON_H
