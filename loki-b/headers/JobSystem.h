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

    class Range {
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

        Range(const Range &other) : isLog(other.isLog), start(other.start), stop(other.stop), n(other.n) {}

        bool next() {
            return (n - 1) > iter++;
        }

        // JvD: also support n=1, a degenerate range with only one value (the start value). Avoid /0.
        double value() const {
            if (isLog)
                return n==1 ? std::pow(10.,start) : std::pow(10., start + iter * (stop - start) / (n - 1));

            return n==1 ? start : start + iter * (stop - start) / (n - 1);
        }

        void reset() {
            iter = 0;
        }
    };

    struct Job {
        std::string name;

        void (WorkingConditions::*callback)(double){nullptr};

        Range range;
    };

    class JobManager {
        std::vector<Job> jobs;
        uint32_t jobIndex{0};

        WorkingConditions *wc;

    public:
        explicit JobManager(WorkingConditions *workingConditions) : wc(workingConditions) {}

        ~JobManager() = default;

        JobManager(const JobManager &other) = delete;

        void addJob(Job &&job) {
            jobs.emplace_back(job);
        }

        void addJob(Job &job) {
            jobs.emplace_back(job);
        }

        void prepareFirstJob() {
            for (Job& job : jobs)
            {
                (wc->*job.callback)(job.range.value());
            }
        }
        bool nextJob() {
            Job &job = jobs[jobIndex];


            if (job.range.next()) {
                (wc->*job.callback)(job.range.value());

                if (jobIndex != jobs.size() - 1)
                    ++jobIndex;

                return true;
            } else {
                if (jobIndex == 0) return false;

                job.range.reset();
                --jobIndex;

                return nextJob();
            }
        }

        std::string getCurrentJobFolder() const {
            std::stringstream ss;

            for (const auto &job : jobs) {
                ss << "_" << job.name << "_" << job.range.value();
            }

            return ss.str();
        }
    };
}


#endif //LOKI_CPP_JOBSYSTEM_H
