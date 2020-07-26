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
#include <cassert>

namespace loki {

    /** This class gives access to a sequence of values and has two public
     *  members: size() returns the number of values, value(size_type)
     *  returns the value for a given index. The number of values is passed
     *  to the constructor and kept in a private member, value() delegates
     *  the retrieval of the value to the protected virtual member get_value.
     *  That member must be implemented by derived classes to do the
     *  actual calculation and retrun the value.
     */
    class Range
    {
    public:
        using size_type = std::size_t;
        /// This constructor records \a size, the number of values in the Range
        Range(size_type size) : m_size(size) {}
        virtual ~Range(){}

        /// returns the number of values of this Range.
        size_type size() const { return m_size; }
        /// Returns the \a ndx'th value. Argument \a ndx must be smaller than size().
        double value(size_type ndx) const
        {
            assert(ndx<size());
            return get_value(ndx);
        }
    protected:
        /// Must be overridden to return the value for \a ndx < size().
        virtual double get_value(size_type ndx) const=0;
    private:
        /// The number of values of this Range
        const size_type m_size;
    };

    /** A Job controls one of the parameters of a parametrized model.
     *  It manages a Range object that describes the value(s) of the parameter
     *  for which the model must be run, In addition it manages a callback function,
     *  which is called by the JobMaanager when this parameter value changes:
     *  it must prepare the model to do a run with the new set of values.
     */
    class Job
    {
    public:
        using callback_type = std::function<void(double)>;
        /** The callback function must accept a double and returns a void.
         */
        Job(const std::string& _name, const callback_type _callback, Range* _range);

        std::string name;
        /** callback is the function that gets called by the JobManager when a
         *  new value in the range is activated.
         */
        const callback_type callback;
        double active_value() const { return range->value(active_ndx); }
        void reset() { active_ndx = 0; }
        bool advance() { return (range->size()-1) > active_ndx++; }
    private:
        std::unique_ptr<const Range> range;
        Range::size_type active_ndx;
    };

    class JobManager
    {
    public:
        JobManager();
        ~JobManager();
        JobManager(const JobManager &other) = delete;

        void addParameter(const std::string& _name, const Job::callback_type _callback, const std::string& range);
        void addParameter(const std::string& _name, const Job::callback_type _callback, const json_type& range);
        void prepareFirstJob();
        bool prepareNextJob();
        std::string getCurrentJobFolder() const;
    private:
        std::vector<std::unique_ptr<Job>> jobs;
        using size_type = std::size_t;
        size_type jobIndex;
    };

} // namespace loki


#endif //LOKI_CPP_JOBSYSTEM_H
