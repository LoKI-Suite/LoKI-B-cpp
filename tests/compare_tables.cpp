#include <iostream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include "LoKI-B/LookupTable.h"

loki::LookupTable extractLUT(std::istream& is)
{
    std::string str((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
    std::stringstream ss{"-----\n" + str + "-----\n"};
    loki::LookupTable lut{loki::LookupTable::create(ss)};
    return lut;
}

double compareTables(const std::string& ifname1, const std::string& ifname2, const std::string& ofname)
{
    std::ifstream ifs1{ifname1};
    if (!ifs1)
    {
        throw std::runtime_error("Could not open file '" + ifname1 + "' for reading.");
    }
    const loki::LookupTable lut_bolsig{extractLUT(ifs1)};
    // now read the loki file
    std::ifstream ifs2{ifname2};
    if (!ifs2)
    {
        throw std::runtime_error("Could not open file '" + ifname2 + "' for reading.");
    }
    std::string header;
    std::getline(ifs2,header);
    std::cout << "#Read loki header: " << header << std::endl;
    const loki::LookupTable lut_loki{extractLUT(ifs2)};

    std::ofstream ofs{ofname};
    if (!ofs)
    {
        throw std::runtime_error("Could not open file '" + ofname + "' for reading.");
    }
    // now make a table based on the Bolsig+ energy points.
    // it prints the energy, the Bolsig+ table value and the loki table value
    // (note that the loki values will feature a table interpolation error when
    // a loki energy point is not also a Bolsig+ energy point).
    double max_rel_diff = 0.0;
    for (loki::LookupTable::Vector::Index i=0; i!= lut_bolsig.size(); ++i)
    {
        const double e = lut_bolsig.x()[i];
        const double v1 = lut_bolsig.interpolate(e);
        const double v2 = lut_loki.interpolate(e);
        ofs << e << '\t' << v1 << '\t' << v2 << std::endl;
        max_rel_diff = std::max(max_rel_diff,std::abs((v1-v2)/v1));
    }
    std::cout << "Combined table written to file '" << ofname << "'." << std::endl;
    return max_rel_diff;
}

const char* usage =
    "Usage: compare_table <file1> <file2> <out-file>.\n"
    "\n"
    "example: ./compare_table eedf_R4.out loki_eedf.out combined.out\n"
    "will construct a lookup table from eedf_R4.out. It writes this lookup table\n"
    "with xtra columns, representing the interpolated values of <file2> and the\n"
    "relative differences.\n";

int main(int argc, const char* argv[])
{
    try {
        if (argc!=4)
        {
            throw std::runtime_error(usage);
        }
        const std::string ifname1(argv[1]);
        const std::string ifname2(argv[2]);
        const std::string ofname(argv[3]);
        const double max_rel_diff = compareTables(ifname1,ifname2,ofname);
        std::cout << "# max_rel_diff: " << max_rel_diff << std::endl;
        return 0;
    }
    catch(std::exception& exc)
    {
        std::cerr << exc.what() << std::endl;
        return 1;
    }
}

