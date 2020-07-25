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

double Range::value() const
{
    if (isLog)
        return n==1 ? std::pow(10.,start) : std::pow(10., start + iter * (stop - start) / (n - 1));
    return n==1 ? start : start + iter * (stop - start) / (n - 1);
}

Job::Job(const std::string& _name, const callback_type _callback, Range* _range)
 : name(_name), callback(_callback), range(_range)
{
}

JobManager::JobManager(WorkingConditions *workingConditions)
 : wc(workingConditions)
{
}

void JobManager::addJob(const std::string& _name, const Job::callback_type _callback, const std::string& range)
{
    jobs.emplace_back(Job{_name, _callback, new Range{range}});
}

void JobManager::addJob(const std::string& _name, const Job::callback_type _callback, const json_type& range)
{
    jobs.emplace_back(Job{_name, _callback, new Range{range}});
}

void JobManager::prepareFirstJob()
{
    for (Job& job : jobs)
    {
        (wc->*job.callback)(job.range->value());
    }
}

bool JobManager::nextJob()
{
    Job &job = jobs[jobIndex];


    if (job.range->next()) {
        (wc->*job.callback)(job.range->value());

        if (jobIndex != jobs.size() - 1)
            ++jobIndex;

        return true;
    } else {
        if (jobIndex == 0) return false;

        job.range->reset();
        --jobIndex;

        return nextJob();
    }
}

std::string JobManager::getCurrentJobFolder() const {
    std::stringstream ss;

    for (const auto &job : jobs) {
        ss << "_" << job.name << "_" << job.range->value();
    }

    return ss.str();
}

} // namespace loki
