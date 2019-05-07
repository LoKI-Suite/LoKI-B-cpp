//
// Created by daan on 6-5-19.
//

#ifndef LOKI_CPP_ALTERNATIVESETUP_H
#define LOKI_CPP_ALTERNATIVESETUP_H

#include <map>
#include <vector>
#include <string>

/*
 * The SetupStruct class provides a way to load a completely arbitrary setup file
 * (that condones to the layout rules).
 */

namespace loki {
    class SetupStruct {
        const std::string inputPath{"Input"};

        uint8_t level;
        std::map<std::string, uint16_t> key;
        std::vector<SetupStruct> vSetupStruct;
        std::vector<std::string> vString;

    public:
        explicit SetupStruct(uint8_t level = 0) : level(level) {}

        ~SetupStruct() = default;

        template <typename T>
        T& add(const std::string &name, const T &value) {}

        // Note that it is unnecessary to make this a template function
        // it might even be preferred to split this function into
        // "getString" and "getSetupStruct" for readability.
        template <typename T>
        T &get(const std::string &name) {}

        SetupStruct &operator[](const std::string &name);

        std::string &operator()(const std::string &name);

        std::vector<std::string> &strings() {return vString;}

        bool parseFile(const std::string &fileName);

        bool parseSection(const std::string &content);
    };
}


#endif //LOKI_CPP_ALTERNATIVESETUP_H
