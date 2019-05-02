#include <iostream>
#include <fstream>
#include <string>
#include <regex>

#define COMMENT_CHAR '%'

void removeComments(std::string &line);

int main() {

    std::ifstream file("../../LoKI/Code/Input/Databases/harmonicFrequencies.txt");

    if (!file.is_open())
    {
        std::cerr << "could not open database file" << std::endl;
        return 1;
    }

    std::regex r(R"((\w*)\s+(.+))");
    std::smatch m;

    std::string line;

    while (getline(file, line))
    {
        removeComments(line);

        if (line.empty()) continue;

        if (std::regex_match(line, m, r))
        {
            std::cout << "Gas: " << m[1] << "\tvalue: " << m[2] << std::endl;
        }
    }

    file.close();

    return 0;
}

void removeComments(std::string &line)
{
    auto it = std::find(line.begin(), line.end(), COMMENT_CHAR);
    line.erase(it, line.end());
}