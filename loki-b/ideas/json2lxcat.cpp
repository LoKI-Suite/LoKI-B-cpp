#include <string>
#include <set>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "nlohmann/json.hpp"

using json_type = nlohmann::json;

/** Make all characters in \a s uppercase; return a reference to the
 *  modified string.
 */
std::string& to_upper_inplace(std::string& s)
{
    for (auto& c : s) c = std::toupper(c);
    return s;
}

/** Make a copy of string \a s, make all characters of the copy
 *  uppercase and return the result (\a s is not modified).
 */
std::string to_upper_copy(const std::string& s)
{
    std::string tmp(s);
    to_upper_inplace(tmp);
    return tmp;
}

/** This function is used only by the implementation of string_to<T>.
 *  It is called after an item of type T is extracted from the stream \a ss
 *  and checks that the stream is not in the fail state and that the eof
 *  is reached. This means that conversion happened, and that no trailing
 *  characters are left. If any of these conditions is not met, a runtime
 *  error is thrown that contains the offending string and the reason
 *  of failure.
 */
void check_conversion(std::istream& ss, const std::string& str, const std::string& type)
{
    if (ss.fail())
    {
        throw std::runtime_error("Error converting '" + str + "' to type '"
            + type + "'.");
    }
    if (!ss.eof())
    {
        std::string s;
        ss >> s;
        if (!s.empty())
        {
            throw std::runtime_error("When converting '" + str + "' to type '"
                + type
                + "' : Trailing characters: '" + s + "'.");
        }
    }
}

/** Convert \a str to a variable of type T and return the result.
 *  A std::runtime_error is thrown if the conversion failed or is
 *  incomplete.
 */
template <class T>
T string_to(const std::string& str)
{
    std::istringstream is(str);
    T val;
    is >> val;
    check_conversion(is, str, typeid(T).name());
    return val;
}


/** Find and return the name of the first state other than "e" that appears
 *  in a stochiometric array \a sa ('one side of a reaction equation').
 *
 *  \todo Check that no second target is present
 */
std::string get_target(const json_type& sa)
{
    for (const auto& s: sa)
    {
        const json_type& state = s.at("state");
        if (state.get<std::string>()!="e")
        {
            return state;
        }
    }
    throw std::runtime_error("Could not identify target state while parsing array '"
            + sa.dump() + "'.");
}

/** In the array \a par_array of strings, locate the first entry that starts with
 *  the string \a pattern and return the result. When \a strip is true, the
 *  leading pattern substring is removed from the result. If no matching entry is
 *  found, a std::runtime_error is thrown.
 */
std::string get_parameter(const json_type& par_array, const std::string& pattern, bool strip)
{
    for (const auto& p : par_array)
    {
        const std::string s(p);
        if (s.substr(0,pattern.size())==pattern)
        {
            return strip ? s.substr(pattern.size()) : s;
        }
    }
    throw std::runtime_error("Parameter starting with '" + pattern + "' not found.");
}

/** Call function get_parameter to get the first entry of the \a par_array that
 *  starts with \a pattern, then convert the remainder of that parameter string
 *  to a double and return the result. If no match is found, or the conversion
 *  fails, a std::runtime_error is thrown.
 *  Example: if the array contains params = [ "A=2","B=3.14" ], the call
 *  get_parameter_value(params,"B=") results in the value 3.14.
 */
double get_parameter_value(const json_type& par_array, const std::string& pattern)
{
    return string_to<double>(get_parameter(par_array,pattern,true));
}

/** The strings that identify process types in LXCat files. The order that
 *  is used here defines a priority order, see the documentation of get_type().
 *
 *  \todo types VIBRATIONAL, ROTATIONAL are not announced in the LXCat headers.
 *  \todo Bolsig+ also appears to support MOMENTUM, the alternative spelling
 *        IONISATION and perhaps more, no idea if we should be prepared for
 *        those as well when parsing LXCat files.
 */
const std::vector<std::string>& process_types()
{
    static const std::vector<std::string> types{
        "ELASTIC",
        "EFFECTIVE",
        "ATTACHMENT",
        "IONIZATION",
        "EXCITATION",
        "VIBRATIONAL",
        "ROTATIONAL"
    };
    return types;
}

/** Check that \a type equals one of the values in process_types().
 *  If this is the case, a reference to \a type is returned, otherwise
 *  a std::runtime_error is thrown.
 */
const std::string& check_type(const std::string& type)
{
    if (std::find(process_types().begin(),process_types().end(),type)==process_types().end())
    {
        throw std::runtime_error("Unknown process type '" + type + "'.");
    }
    return type;
}

/** Select a type from the \a type_tags array and use that in the LXCat file,
 *  do translation if necessary and write the result in uppercase.
 *
 *  Note that JSON files allow for more than one tag, e.g. [ "Electronic",
 *  "Vibrational" ]. This function will then prefer the type that appears
 *  first in the list of allowed strings that is returned by process_types().
 *  As a special rule, the string Electronic, which identifies an electronic
 *  excitation, is translated into EXCITATION.
 *
 *  \todo establish and document a clear set of rules.
 */
std::string get_type(const json_type& type_tags)
{
    std::set<std::string> tags;
    for (const auto& t : type_tags)
    {
        tags.insert(t=="Electronic" ? "EXCITATION" : to_upper_copy(t));
    }
    if (tags.size()==0)
    {
        throw std::runtime_error("At least one process tag must be specified.");
    }
    else if (tags.size()==1)
    {
        return check_type(*tags.begin());
    }
    else
    {
#if 0
        os << "COMMENT: type tags in the JSON file: ";
        for (const auto& t : tags)
        {
            std::cout << ' ' << t;
        }
        std::cout << std::endl;
#endif
        /// \todo What is the 'preferred' process type if more are given?
        for (const auto& t : process_types())
        {
            if (tags.find(t)!=tags.end())
            {
                return t;
            }
        }
        throw std::runtime_error("Unknoqn process type(s).");
    }
}

/** Write a cross section table that is defined by the JSON object \a table
 *  to the stream \a os. The two axis lables and units are obtained from
 *  arrays "labels" and "units" in \a table, each must contain two strings.
 *  This information is used to produce the "COLUMNS:" line in the LXCat file.
 *  the xy-data are obtained from the array "data" in \a table, which must
 *  contain arrays of two doubles that each represent a data point.
 */
void emit_table(const json_type& table, std::ostream& os)
{
    const json_type& labels = table.at("labels");
    if (labels.size()!=2)
    {
        throw std::runtime_error("Expected two labels.");
    }
    const json_type& units = table.at("units");
    if (units.size()!=2)
    {
        throw std::runtime_error("Expected two units.");
    }
    os << "COLUMNS: "
        << labels[0].get<std::string>() << " (" << units[0].get<std::string>() << ")"
            << " | "
        << labels[1].get<std::string>() << " (" << units[1].get<std::string>() << ")"
            << '\n';
    os << "------" << '\n';
    /// \todo Check that the table is not empty? Do we need one or two entries at least?
    for (const auto& entry : table.at("data"))
    {
        if (entry.size()!=2)
        {
            throw std::runtime_error("Expected two values in a table entry.");
        }
        os << ' ' << entry[0].get<double>() << '\t' << entry[1].get<double>() << '\n';
    }
    os << "------" << '\n';
}

/** Write one side of a reaction equation, as defined by the object \a sa,
 *  to the stream \a os. The object \a sa must be an array of objects, each
 *  containing "state" (string) and "count" (int) values, representing the
 *  name and occurrence of a particular state in the format. As an example,
 *  the state array
 *
 *    [ { "state": "Ar+", "count": 1 }, { "state": "e", "count": 2 } ]
 *
 *  would result in the output "Ar+ + 2 e". Note that coefficient values
 *  equal to unity are not written. When \a skip_e is true, entries
 *  with "state" equal to "e" (representing an electron) are skipped.
 */
void emit_stoich_array(const json_type& sa, std::ostream& os, bool skip_e)
{
    bool first_item=true;
    for (const auto& s: sa)
    {
        const std::string sname = s.at("state").get<std::string>();
        if (skip_e && sname=="e")
        {
            continue;
        }
        if (first_item)
        {
            // first item. No special handling needed.
            first_item = false;
        }
        else
        {
            os << " + ";
        }
        unsigned c = s.contains("count") ? s.at("count").get<unsigned>() : 1;
        if (c!=1)
        {
            os << c << ' ';
        }
        os << sname;
    }
}

/** write the reaction equation that is defined by \a fmt to \a os.
 *  This produces the output "<LHS> <SEP> <RHS>", where the sides
 *  LHS and RHS are produced by calling emit_stoich_array on the
 *  objects "lhs" and "rhs" of \a fmt, respectively. The separator
 *  is "<->" when "reversible" has value true, "->" otherwise.
 *  The argument \a skip_e is passed on to the calls to emit_stoich_array,
 *  when equal to true the electron entries are skipped for both the
 *  left and right-hand side of the reaction equation.
 *  The function returns the value of "reversible".
 */
bool emit_format(const json_type& fmt, std::ostream& os, bool skip_e)
{
    emit_stoich_array(fmt.at("lhs"),os,skip_e);
    os << ' ';
    const bool reversible = fmt.at("reversible");
    if (reversible)
        os << "<->";
    else
        os << "->";
    os << ' ';
    emit_stoich_array(fmt.at("rhs"),os,skip_e);
    os << '\n';
    return reversible;
}

/** Write information that is particular for a process of \a type ELASTIC
 *  or EFFECTIVE, defined by the object \a pnode, to \a os.
 *
 *  For a process "e + X -> e + X", the code writes the collision target
 *  (here X), followed by a line contain the mass ratio m_e/m_x. The latter
 *  is extracted from the parameter in the "parameters" string list that is
 *  of the form "m/M = <value>".
 */
void emit_elastic_effective(const json_type& pnode, const std::string& type, std::ostream& os)
{
    const json_type& rnode = pnode.at("reaction");
    os << get_target(rnode.at("lhs")) << '\n';
    /** \todo Maybe let the parameter array in the JSON be an array of objects so we
     *        can look for attribute "m/M".
     */
    os << ' ' << get_parameter_value(pnode.at("parameters"),"m/M = ") << '\n';
}

/** Write information that is particular for inelastic processes (electronic,
 *  vibrational and rotational excitation).
 *
 *  Firstly, the reaction format is written, but without the electrons. For an excitation
 *  or ionization process this results in something like "X -> X*" or "X -> X+". The
 *  subsequent line contains the energy threshold of the process in eV (just the number,
 *  without the units). When the process is marked to be reversible in the JSON file
 *  '<->' is written instead of '->'. If a reversible process is an electronic excitation,
 *  the ratio of the electronic statistical weights is written on the line following
 *  that containing the threshold. The value is obtained from the paramater string that
 *  starts with "g1/g0 = ".
 *
 *  \todo Other processes, like attachment, perhaps need different handling? Check.
 *  \todo Are parameters available for every inelastic process? Should we implement
 *        the handling of parameters differently for each process type?
 */
void emit_inelastic(const json_type& pnode, const std::string& type, std::ostream& os)
{
    const bool skip_e = true;
    bool reversible = emit_format(pnode.at("reaction"),os, skip_e);
    os << ' ' << pnode.at("threshold") << '\n' ;
    /** \todo When is the g1/g0 required? In the N2 files there are processes
     *        that have <->, also have the g-ration in the parameter list, but
     *        not on 'line 3' (after the energy threshold). Is this needed
     *        only for electronic excitation?
     */
    if (reversible && type=="EXCITATION")
    {
        os << ' ' << get_parameter_value(pnode.at("parameters"),"g1/g0 = ") << '\n';
    }
}

/** \todo document me.
 */
void emit_process(const json_type& pnode, std::ostream& os)
{
    const json_type& rnode = pnode.at("reaction");
    const std::string type = get_type(rnode.at("type_tags"));

    // 1. Emit the process type

    os << type << '\n';

    // 2. Emit process-type-specific bits

    if (type=="ELASTIC" || type=="EFFECTIVE")
    {
        /// \todo double-check that EFFECTIVE is handled the same way as ELASTIC
        emit_elastic_effective(pnode,type,os);
    }
    else
    {
        emit_inelastic(pnode,type,os);
    }

    // 3. Emit common information:
    //    - parameters
    //    - references
    //    - the cross section table

    /** Generic handling of parameters. Produce a line that starts
     *  with the word "PARAM. ", followed by a comma-separated list
     *  of the parameter string values. Depending on the process type,
     *  the parameter info also ends up elswhere in the LXCat output
     *  (m/M for elastic processes, g1/g0 for reversible electron
     *  excitation, etc.).
     *
     *  \todo Find out and document which information must be declared
     *        for each process type, and make it so.
     *  \todo Some databases have systematic ways of providing additional
     *        information in the form of comments. As an example, the IST
     *        databases have the full reaction with detailed state info
     *        in addition to the simplified format. We may need a program
     *        option like "--style=IST" to announce that an IST file is
     *        being processed, so we can add thos extras as well. This may
     *        require us to read the state list as well, because that can
     *        give us things like v,J quantum numbers for a state that we
     *        know only by name otherwise.
     */
    if (pnode.contains("parameters"))
    {
        os << "PARAM.: ";
        bool first=true;
        for (const auto& p : pnode.at("parameters"))
        {
            if (first)
            {
                first = false;
            }
            else
            {
                os << ", ";
            }
            os << p.get<std::string>();
        }
        os << '\n';
    }

    /** \todo Are there any limits on the line size? Should lines be,
     *   clipped to 80 characters, for example?
     */
    for (const auto& r : pnode.at("reference"))
    {
        os << "COMMENT: " << r.get<std::string>() << '\n';
    }
    emit_table(pnode,os);
    os << '\n';
}

/** No string preference what this should be. 80 is 'classic',
 *  but somewhat small on terminals that are less than 25 years old.
 */
const std::string::size_type line_length = 100;

void format_meta(std::ostream& os, const std::string& prefix, std::string str)
{
    bool first=true;
    const std::string::size_type part_length = line_length - prefix.size();
    while (!str.empty())
    {
        const std::string::size_type nl_pos = str.find("\n");
        if (first)
        {
            os << prefix;
            first = false;
        }
        else
        {
            os << std::string(prefix.size(),' ');
        }
        if (str.size()<=part_length && nl_pos==std::string::npos)
        {
            os << str;
            str = std::string{};
        }
        else
        {
            const std::string::size_type pos = std::min(nl_pos,str.find_first_of(" \t",70));
            os << str.substr(0,pos);
            str = str.substr(pos+1);
        }
        os << std::endl;
    }
}

/** Write the JSON (v6) document \a src to \a os in LXCat format.
 *
 *  \todo The general LXCat header is not yet written.
 *  \todo That header is incomplete, it seems. In the IST LXCat file
 *        the types VIBRATIONAL and ROTATIONAL are not mentioned. What
 *        is supported, and what are the process-specific expectations.
 *  \todo The database-specific header is not written. As an example,
 *        take the file Nitrogen/N2_LXCat.txt. There we see that this
 *        comes in two parts: at the top we have "RECOMMENDED REFERENCE
 *        FORMAT", followed by an IST-block. Then we have, after the
 *        line with xxxxxxx the DATABASE, PERMLINK, DESCRIPTION, CONTACT
 *        and HOW TO REFERENCE lines (are these really standard?), then
 *        there is the data set specific (?) info in the form of
 *        COMMENT lines *after the *****. This is different in the
 *        vib and rot files, so we cannot speak of a common "IST"
 *        block.
 *  \todo In an LXCat file, you would find, for example, N2 as name
 *        for the ground state of nitrogen, in the JSON files we now
 *        have N2{X^1S_g^+}. is that a problem? It definitely is not
 *        compatiable. Should we (optionally) translate the latter in
 *        the former?
 *  \todo In LXCat you typically see Gas(State), in JSON we have Gas{State}.
 *        Existing users may rely on parentheses being used to announce
 *        a state.
 */
void json2lxcat(const json_type& src, std::ostream& os)
{
    unsigned ndx=0;
    for (const auto& sit: src.at("sets"))
    {
        os << std::string(line_length,'*') << '\n' << '\n';
        format_meta(os,"COMMENT: ",sit.at("description"));
        os << '\n';
        os << std::string(line_length,'*') << '\n' << '\n';
        for (const auto& pcnf: sit.at("processes"))
        {
            try {
                emit_process(pcnf,os);
                ++ndx;
            }
            catch (std::exception& exc)
            {
                throw std::runtime_error("While parsing process #"
    		+ std::to_string(ndx) + ": " + std::string(exc.what()));
            }
        }
    }
}

int main()
{
    try {
        const std::string str(
            std::istreambuf_iterator<char>{std::cin},
            std::istreambuf_iterator<char>{});
        const json_type src(json_type::parse(str));
        json2lxcat(src,std::cout);
        return 0;
    }
    catch (std::exception& exc)
    {
        std::cerr << exc.what() << std::endl;
        return 1;
    }
}
