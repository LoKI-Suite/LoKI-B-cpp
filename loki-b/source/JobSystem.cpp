//
// Created by daan on 08-07-2019.
//

#include "JobSystem.h"
#include "Parse.h"

namespace loki {

namespace impl {

    static Range getRange(const std::string &rangeString, bool &success)
    {
        static const std::regex r(
            R"(\s*((?:logspace\()|(?:linspace\())\s*(-?\d+\.?\d*)\s*,\s*(-?\d+\.?\d*)\s*,\s*(\d+\.?\d*))");
        std::smatch m;

        // No checking since this has already been performed once this function is called.
        if (!std::regex_search(rangeString, m, r))
        {
            success = false;
            return Range(0., 0., 0, false);
        }

        std::stringstream ss;

        const std::string function = m.str(1);

        double start, stop;
        uint32_t steps;

        ss << m[2];
        ss >> start;
        ss.clear();
        ss << m[3];
        ss >> stop;
        ss.clear();
        ss << m[4];
        ss >> steps;

        success = true;

        return Range(start, stop, steps, function[1] == 'o');
    }

} // namespace impl

Range::Range(const std::string& str)
{
    if (Parse::isNumerical(str))
    {
        double value;
        bool success = Parse::getValue(str,value);
        if (!success) {
            throw std::runtime_error("Invalid value/range specification '" + str + "'.");
        }
        *this = Range{value};
    }
    else
    {
        bool success = false;
        *this = impl::getRange(str, success);
        if (!success) {
            throw std::runtime_error("Invalid value/range specification '" + str + "'.");
        }
    }
}

Range::Range(const json_type& cnf)
{
    if (cnf.type()==json_type::value_t::string)
    {
        bool success = true;
        *this = impl::getRange(cnf, success);
        if (!success) {
            throw std::runtime_error("Invalid value/range specification '" + cnf.dump(2) + "'.");
        }
    }
    else
    {
        try {
            *this = Range(cnf.get<double>());
        }
        catch(std::exception& exc)
        {
            throw std::runtime_error("Invalid value/range specification '" + cnf.dump(2) + "': "
                    + std::string(exc.what()));
        }
    }
}

} // namespace loki
