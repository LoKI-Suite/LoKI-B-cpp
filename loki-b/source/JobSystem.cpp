//
// Created by daan on 08-07-2019.
//

#include "JobSystem.h"
#include "Parse.h"

namespace loki {

/** A degenerate Range implementation that describes a single value.
 */
class RangeSingleValue : public Range
{
public:
    RangeSingleValue(double value)
        : Range(1), m_value(value)
    {
    }
protected:
    virtual double get_value(size_type ndx) const
    {
        return m_value;
    }
private:
    double m_value;
};

/** A range that represents a 'linspan', a series of equidistant
 *  values in a given domain.
 */
class RangeLinSpan : public Range
{
public:
    RangeLinSpan(double start, double stop, size_type steps)
        : Range(steps), start(start), stop(stop)
    {
        if (size()<2)
        {
            throw std::runtime_error("Range(linspan): at least two points required.");
        }
    }
protected:
    virtual double get_value(size_type ndx) const
    {
        return start + ndx * (stop - start) / (size() - 1);
    }
private:
    double start;
    double stop;
};

/** A range that represents a 'logspan', the values are obtained as
 *  10^v, where v is a series of equidistant values in a given
 *  domain.
 */
class RangeLogSpan : public Range
{
public:
    RangeLogSpan(double log_start, double log_stop, size_type steps)
        : Range(steps), log_start(log_start), log_stop(log_stop)
    {
        if (size()<2)
        {
            throw std::runtime_error("Range(logspan): at least two points required.");
        }
    }
protected:
    virtual double get_value(size_type ndx) const
    {
        const double log_value = log_start + ndx * (log_stop - log_start) / (size() - 1);
        return std::pow(10.,log_value);
    }
private:
    double log_start;
    double log_stop;
};

static Range* createLinLogRange(const std::string &rangeString)
{
    try {
        static const std::regex r(
            R"(\s*((?:logspace\()|(?:linspace\())\s*(-?\d+\.?\d*)\s*,\s*(-?\d+\.?\d*)\s*,\s*(\d+\.?\d*))");
        std::smatch m;

        if (!std::regex_search(rangeString, m, r))
        {
            throw std::runtime_error("Parse error.");
        }

        std::stringstream ss;

        const std::string function = m.str(1);

        double start, stop;
        Range::size_type steps;

        ss << m[2];
        ss >> start;
        ss.clear();
        ss << m[3];
        ss >> stop;
        ss.clear();
        ss << m[4];
        ss >> steps;

        if (function=="linspan")
            return new RangeLinSpan(start, stop, steps);
        else if (function=="logspan")
            return new RangeLogSpan(start, stop, steps);
        else
            throw std::runtime_error("Unknown function '" + function + "'.");
    }
    catch (std::exception& exc)
    {
        throw std::runtime_error("Error creating range from string '" + rangeString + ":\n"
            + std::string(exc.what()));
    }
}

Range* createRange(const std::string& str)
{
    if (Parse::isNumerical(str))
    {
        double value;
        bool success = Parse::getValue(str,value);
        if (!success) {
            throw std::runtime_error("Invalid value/range specification '" + str + "'.");
        }
        return new RangeSingleValue{value};
    }
    else
    {
        return createLinLogRange(str);
    }
}

Range* createRange(const json_type& cnf)
{
    if (cnf.type()==json_type::value_t::string)
    {
        return createLinLogRange(cnf);
    }
    else
    {
        return new RangeSingleValue(cnf.get<double>());
    }
}

Job::Job(const std::string& _name, const callback_type _callback, Range* _range)
 : name(_name), callback(_callback), range(_range), active_ndx(0)
{
}

JobManager::JobManager()
 : jobIndex(0)
{
}

JobManager::~JobManager()
{
}

void JobManager::addParameter(const std::string& _name, const Job::callback_type _callback, const std::string& range)
{
    jobs.emplace_back(new Job{_name, _callback, createRange(range)});
}

void JobManager::addParameter(const std::string& _name, const Job::callback_type _callback, const json_type& range)
{
    jobs.emplace_back(new Job{_name, _callback, createRange(range)});
}

void JobManager::prepareFirstJob()
{
    for (const auto& parameter : jobs)
    {
        (parameter->callback)(parameter->active_value());
    }
}

bool JobManager::prepareNextJob()
{
    Job &parameter = *jobs[jobIndex];


    if (parameter.advance()) {
        (parameter.callback)(parameter.active_value());

        if (jobIndex != jobs.size() - 1)
            ++jobIndex;

        return true;
    } else {
        if (jobIndex == 0) return false;

        parameter.reset();
        --jobIndex;

        return prepareNextJob();
    }
}

std::string JobManager::getCurrentJobFolder() const {
    std::stringstream ss;

    for (const auto &parameter : jobs) {
        ss << "_" << parameter->name << "_" << parameter->active_value();
    }

    return ss.str();
}

} // namespace loki
