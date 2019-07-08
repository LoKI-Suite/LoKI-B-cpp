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

namespace loki {
    struct WorkingConditions;

    class Range {
        const bool isLog;
        const double start, stop;

        const uint32_t n;
        uint32_t iter{0};

    public:
        Range(double start, double stop, uint32_t nSteps, bool isLog) : start(start), stop(stop), n(nSteps),
                                                                        isLog(isLog) {}

        Range(const Range &other) : start(other.start), stop(other.stop), n(other.n), isLog(other.isLog) {}

        bool next() {
            return (n - 1) > iter++;
        }

        double value() const {
            if (isLog)
                return std::pow(10., start + iter * (stop - start) / (n - 1));

            return start + iter * (stop - start) / (n - 1);
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
