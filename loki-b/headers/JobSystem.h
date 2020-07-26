//
// Created by daan on 08-07-2019.
//

#ifndef LOKI_CPP_JOBSYSTEM_H
#define LOKI_CPP_JOBSYSTEM_H

#include <cstdint>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <functional>
#include "json.h"

namespace loki {

    class Range
    {
        /** \todo make isLog, start, stop and n constant after the string and json_type ctors
         *  have been updated to allow that.
         */
        bool isLog;
        double start, stop;

        uint32_t n;
        uint32_t iter{0};
    public:
        // expects a legacy Range string: a numerical value or a linspan or logspan description
        Range(const std::string& str);
        Range(const json_type& cnf);
        Range(double value) : isLog(false), start(value), stop(value), n(1) {}
        Range(double start, double stop, uint32_t nSteps, bool isLog) : isLog(isLog), start(start), stop(stop), n(nSteps) {}

        void reset() { iter = 0; }
        bool next() { return (n-1) > iter++; }
        double value() const;
    };

    /** A Job controls one of the parameters of a parametrized model.
     *  It manages a Range object that describes the value(s) of the parameter
     *  for which the model must be run, In addition it manages a callback function,
     *  which is called by the JobMaanager when this parameter value changes:
     *  it must prepare the model to do a run with the new set of values.
     */
    struct Job
    {
        using callback_type = std::function<void(double)>;
        /** The callback function must accept a double and returns a void.
         */
        Job(const std::string& _name, const callback_type _callback, Range* _range);

        std::string name;
        /** callback is the function that gets called by the JobManager when a
         *  new value in the range is activated.
         */
        const callback_type callback;
        std::unique_ptr<Range> range;
    };

    class JobManager
    {
    public:
        JobManager();
        ~JobManager() = default;
        JobManager(const JobManager &other) = delete;

        void addJob(const std::string& _name, const Job::callback_type _callback, const std::string& range);
        void addJob(const std::string& _name, const Job::callback_type _callback, const json_type& range);
        void prepareFirstJob();
        bool nextJob();
        std::string getCurrentJobFolder() const;
    private:
        std::vector<Job> jobs;
        uint32_t jobIndex{0};
    };

} // namespace loki


#endif //LOKI_CPP_JOBSYSTEM_H
