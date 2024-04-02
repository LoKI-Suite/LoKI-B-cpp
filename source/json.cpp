#include "LoKI-B/json.h"
#include "LoKI-B/Log.h"
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace loki
{

json_type read_json_from_stream(std::istream &is)
{
    std::string str((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
    return json_type::parse(str);
}

json_type read_json_from_file(const std::filesystem::path &fname)
{
    std::ifstream is(fname);
    if (!is)
    {
        Log<Message>::Error("Could not open file '", fname.generic_string(), "' for reading.");
    }
    try
    {
        return read_json_from_stream(is);
    }
    catch (std::exception &exc)
    {
        throw std::runtime_error("Error reading file '" + fname.generic_string() + "':\n" + std::string(exc.what()));
    }
}

} // namespace loki
