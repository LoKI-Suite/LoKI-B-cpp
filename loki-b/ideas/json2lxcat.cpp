#include <string>
#include <set>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "nlohmann/json.hpp"

using json_type = nlohmann::json;

std::string& to_upper_inplace(std::string& s)
{
    for (auto& c : s) c = std::toupper(c);
    return s;
}

std::string to_upper_copy(const std::string& s)
{
    std::string tmp(s);
    to_upper_inplace(tmp);
    return tmp;
}

/// \todo Check that no second target is present
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

template <class T>
T string_to(const std::string& str)
{
    std::istringstream is(str);
    T val;
    is >> val;
    check_conversion(is, str, typeid(T).name());
    return val;
}

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

double get_parameter_value(const json_type& par_array, const std::string& pattern)
{
    return string_to<double>(get_parameter(par_array,pattern,true));
}

const std::vector<std::string>& process_types()
{
    /// \todo types VIBRATIONAL, ROTATIONAL are not announced in the header.
    static const std::vector<std::string> types{
        "ELASTIC",
        "EFFECTIVE",
        "ATTACHMENT"
        "IONIZATION",
        "EXCITATION",
        "VIBRATIONAL",
        "ROTATIONAL",
    };
    return types;
}

std::string get_type(const json_type& type_tags)
{
    std::set<std::string> tags;
    for (const auto& t : type_tags)
    {
        tags.insert(to_upper_copy(t));
    }
    if (tags.size()==1)
    {
        return *tags.begin();
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

/** Emit the cross section table.
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

/* 
 * EFFECTIVE
 * N2
 *  1.959210e-5
 * SPECIES: e / N2
 * PROCESS: E + N2 -> E + N2, Effective
 * PARAM.:  m/M = 0.0000195921, complete set
 * COMMENT: [e + N2(X) -> e + N2(X), Effective] Pitchford L C and Phelps A V 1982 Bull. Am. Phys.
 * COMMENT: Soc. 27 109 Tachibana K and Phelps A V 1979 JCP 71 3544.
 * UPDATED: 2017-11-14 10:24:49
 * COLUMNS: Energy (eV) | Cross section (m2)
 * -----------------------------
 *  0.000000e+0    1.100000e-20
 *     ...
 * -----------------------------
 *
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

void emit_elastic_effective(const json_type& pnode, const std::string& type, std::ostream& os)
{
    const json_type& rnode = pnode.at("reaction");
    os << get_target(rnode.at("lhs")) << '\n';
    /** \todo Maybe let the parameter array in the JSON be an array of objects so we
     *        can look for attribute "m/M".
     */
    os << ' ' << get_parameter_value(pnode.at("parameters"),"m/M = ") << '\n';
}

void emit_inelastic(const json_type& pnode, const std::string& type, std::ostream& os)
{
    const bool skip_e = true;
    bool reversible = emit_format(pnode.at("reaction"),os, skip_e);
    os << ' ' << pnode.at("threshold") << '\n' ;
    /** \todo When is the g1/g0 required? In the N2 files there are processes
     *        that have <->, also have the g-ration in the parameter list, but
     *        not on 'line 3' (after the energy threshold). Is this needed
     *        only for electronic exceitation?
     */
    if (reversible && type=="EXCITATION")
    {
        os << ' ' << get_parameter_value(pnode.at("parameters"),"g1/g0 = ") << '\n';
    }
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
    //    - references
    //    - the cross section table

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

void json2lxcat(const json_type& src, std::ostream& os)
{
    for (const auto& pcnf: src.at("processes"))
    {
        emit_process(pcnf,os);
    }
}

int main()
{
    try {
        const std::string str(
            std::istreambuf_iterator<char>(std::cin),
            std::istreambuf_iterator<char>());
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
