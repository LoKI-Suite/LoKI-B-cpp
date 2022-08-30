#include "LoKI-B/LookupTable.h"
#include <stdexcept>

namespace loki {

LookupTable::LookupTable(const Vector& x, const Vector& y)
    : m_x(x), m_y(y)
{
    if (m_x.size()!=m_y.size())
    {
        throw std::runtime_error("Dimensions of abscissa and ordinate do not match.");
    }
    if (m_x.size()<2)
    {
        throw std::runtime_error("At least two points are required.");
    }
    for (Index i=1; i!=m_x.size(); ++i)
    {
        if (m_x[i]<=m_x[i-1])
        {
            throw std::runtime_error("Abscissa values must appear in ascending order.");
        }
    }
}

LookupTable LookupTable::create(const json_type& data)
{
    using PairVector = std::vector<std::pair<double, double>>;
    const PairVector pairs(data.at("data").get<PairVector>());
    Vector x,y;
    x.resize(pairs.size());
    y.resize(pairs.size());
    for (PairVector::size_type i = 0; i != pairs.size(); ++i)
    {
        x[i] = pairs[i].first;
        y[i] = pairs[i].second;
    }
    return LookupTable(x,y);
}

LookupTable LookupTable::create(std::istream& is)
{
    std::vector<double> x, y;
    std::string line;

    while (std::getline(is, line))
    {
        if (line.substr(0, 2) == "--")
        {
            break;
        }
    }
    while (std::getline(is, line))
    {
        if (line.size() && line[0]=='#')
        {
            continue;
        }
        std::stringstream ss(line);
        double energy, value;
        ss >> energy >> value;
        if (!ss)
        {
            // This line did not start with two valid numbers.
            // Check that we indeed reached the end-of-table
            // marker "--" (there may be more "-").
            if (line.substr(0, 2) != "--")
            {
                throw std::runtime_error("Bad data line '" + line + "': expected numbers or '--'");
            }
            // OK, -- was found. Stop reading.
            break;
        }
        x.emplace_back(energy);
        y.emplace_back(value);
    }
    return LookupTable(
            Vector::Map(x.data(), x.size()),
            Vector::Map(y.data(), y.size())
    );
}

} // namespace loki
