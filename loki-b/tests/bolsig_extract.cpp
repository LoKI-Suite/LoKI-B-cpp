#include <iostream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <fstream>

void extractTable(std::istream& is, const std::string& run_id, const std::string& key, std::ostream& os)
{
    std::string line;
    // read until we find the line that starts with run_id.
    // when we find it, print that line with "# Run " in front
    while (!is.eof())
    {
        std::getline(is,line);
        line.erase(line.find_last_not_of(" \n\r\t")+1);
        std::stringstream ss{line};
        std::string id;
        ss >> id;
        if (ss && id==run_id && ss.eof())
        {
            // we just read the case ID, like 'R2'
            os << "# Run '" << run_id << "'" << std::endl;
            break;
        }
    }
    // read until we find a line that contains the key.
    // when we find it, print that line with "# " in front,
    // then read the line that follows (it should contain -----),
    // ignore that line.
    while (!is.eof())
    {
        std::getline(is,line);
        if (line.find(key)!=std::string::npos)
        {
            line.erase(line.find_last_not_of(" \n\r\t")+1);
            os << "# " << line << std::endl;
            // read the line ----- that follows
            std::getline(is,line);
            break;
        }
    }
    // copy lines from the input until an empty line is found
    while (!is.eof())
    {
        std::getline(is,line);
        line.erase(line.find_last_not_of(" \n\r\t")+1);
        if (line.empty())
        {
            break;
        }
        os << line << std::endl;
    }
}

void extractTable(const std::string& ifname, const std::string& run_id, const std::string& key, const std::string& ofname)
{
        std::ifstream ifs{ifname};
        if (!ifs)
        {
            throw std::runtime_error("Could not open file '" + ifname + "' for reading.");
        }
        std::ofstream ofs{ofname};
        if (!ofs)
        {
            throw std::runtime_error("Could not open file '" + ofname + "' for writing.");
        }
        extractTable(ifs,run_id,key,ofs);
}

const char* usage =
    "Usage: bolsig_extract <file_in> <run_id> <key> <file_out>.\n"
    "\n"
    "example: ./bolsig_extract ../loki-b/tests/bolsig_extract.in R4 EEDF eedf_R4.out\n"
    "will extract the EEDF from run R4 from bolsig_extract.in and write it to eedf_R4.out.\n"
    "\n"
    "NOTE: Bolsig+ output must have been generated with the run by run option.\n";

int main(int argc, const char* argv[])
{
    try {
        if (argc!=5)
        {
            throw std::runtime_error(usage);
        }
        const std::string ifname(argv[1]);
        const std::string id(argv[2]);
        const std::string key(argv[3]);
        const std::string ofname(argv[4]);
        extractTable(ifname,id,key,ofname);
        //extractLUT(std::cin,id,key);
        return 0;
    }
    catch(std::exception& exc)
    {
        std::cerr << exc.what() << std::endl;
        return 1;
    }
}
