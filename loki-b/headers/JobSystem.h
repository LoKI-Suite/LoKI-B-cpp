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

    /** A Range class gives access to a sequence of values and has two public
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

        /** Create a new Range object from string \a str. The caller must
         *  assume ownership of the pointer that is returned by this function.
         *
         *  The argument can describe a single value, or a linear of logarithmic
         *  value range. Some sample input and the values it will produce are:
         *  \verbatim
              "42.0"             # 42.0
              "linspan(0,20,3)"  # 0, 10, 20
              "logspan(2,4,3)"   # 1e2, 1e3, 1e4 \endverbatim
         *
         *  \sa RangeSingleValue
         *  \sa RangeLinSpace
         *  \sa RangeLogSpace
         */
        static Range* create(const std::string& str);
        static Range* create(const json_type& cnf);
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
        /// This constructor records \a size, the number of values in the Range
        Range(size_type size) : m_size(size) {}
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
        Job(const std::string& _name, const callback_type _callback, const Range* _range);

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

    /** Class JobManager makanages a collection of parameter definitions and
     *  makes it easy to run a simulation for each combination of parameter values.
     *
     *  A parameter can be declared by calling member addParameter, passing the name
     *  of the parameter, a pointer to a callback function, which is called when
     *  the JobManager activates a new value for this parameter, and a pointer to
     *  a Range object, which describes the vaues of this parameter. This class
     *  takes ownership of the Range pointer and will delete it at the end of its
     *  lifetime.
     */
    class JobManager
    {
    public:
        JobManager();
        ~JobManager();
        JobManager(const JobManager &other) = delete;

        void addParameter(const std::string& _name, const Job::callback_type _callback, const Range* range);
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
