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
#include "json.h"

namespace loki {

    class WorkingConditions;

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
        Range(const Range &other) : isLog(other.isLog), start(other.start), stop(other.stop), n(other.n) {}
        // expects a legacy Range string: a numerical value or a linspan or logspan description
        Range(const std::string& str);
        Range(const json_type& cnf);
        Range(double value) : isLog(false), start(value), stop(value), n(1) {}
        Range(double start, double stop, uint32_t nSteps, bool isLog) : isLog(isLog), start(start), stop(stop), n(nSteps) {}

        void reset() { iter = 0; }
        bool next() { return (n-1) > iter++; }
        double value() const;
    };

    struct Job
    {
        /** Type callback_type is a pointer to a member of WorkingConditions
         *  that accepts a double and returns void.
         */
        using callback_type = void (WorkingConditions::*)(double);
        Job(const std::string& _name, const callback_type _callback, const Range& _range);

        std::string name;
        /** callback is the function that gets called by the JobManager when a
         *  new value in the range is activated.
         */
        const callback_type callback;
        Range range;
    };

    class JobManager
    {
    public:
        explicit JobManager(WorkingConditions *workingConditions);
        ~JobManager() = default;
        JobManager(const JobManager &other) = delete;

        void addJob(Job &&job);
        void addJob(Job &job);
        void prepareFirstJob();
        bool nextJob();
        std::string getCurrentJobFolder() const;
    private:
        std::vector<Job> jobs;
        uint32_t jobIndex{0};
        WorkingConditions *wc;
    };

} // namespace loki


#endif //LOKI_CPP_JOBSYSTEM_H
