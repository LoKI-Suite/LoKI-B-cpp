#include <iostream>
#include <string>
#include <assert.h>
#include <Parse.h>

#define test(x) x(); notify("Test \"" #x "\" was successful.")

using std::cout;
using std::endl;

using namespace loki;

void parse_tests_phase_one();
void parse_tests_phase_two();

void notify(const char *message) {
    std::cout << message << std::endl;
}

int main(int argc, char **argv) {
    notify("Testing: Parsing Phase One");
    parse_tests_phase_one();
    notify("Testing: Parsing Phase Two");
    parse_tests_phase_two();

    notify("All tests complete");
}

/*** PARSING TESTS ***/

/** PHASE ONE **/

void parse_bool() {
    bool value = false;
    std::string field = "test_value: true";
//    std::string field_other = "test_value: TRUE";

    assert(Parse::setField(field, "test_value", value));
    assert(value);

//    assert(Parse::setField(field_other, "test_value", value));
//    assert(value);
}

void parse_double() {
    double value = 0;
    std::string field = "test_value: \t  \t\t  42   \n\n ";


    bool result = Parse::setField(field, "test_value", value);

    assert(result);
    assert(value == 42.);
}

void parse_string() {
    std::string value;
    std::string field = "    electricQuadrupoleMoment: Databases/quadrupoleMoment.txt";
    std::string field_other = "    electricQuadrupoleMoment: \t \tDatabases/quadrupoleMoment.txt  \n\n\n\n  ";
    std::string field_range = " electronTemperature: logspace(11.4574, 23.5, 8) ";
    std::string section = "  gasProperties:                       \n"
                          "    mass: Databases/masses.txt\n"
                          "    fraction:\n"
                          "      - N2 = 1\n"
                          "    harmonicFrequency: Databases/harmonicFrequencies.txt\n"
                          "    anharmonicFrequency: Databases/anharmonicFrequencies.txt\n"
                          "    rotationalConstant: Databases/rotationalConstants.txt\n"
                          "    electricQuadrupoleMoment: Databases/quadrupoleMoment.txt\n"
                          "    OPBParameter: Databases/OPBParameter.txt\n"
                          "  stateProperties:                     \n"
                          "    energy:";

    assert(Parse::setField(field, "electricQuadrupoleMoment", value));
    assert(value == "Databases/quadrupoleMoment.txt");

    assert(Parse::setField(field_other, "electricQuadrupoleMoment", value));
    assert(value == "Databases/quadrupoleMoment.txt");

    assert(Parse::setField(field_range, "electronTemperature", value));
    assert(value == "logspace(11.4574, 23.5, 8)");

    assert(Parse::setField(section, "electricQuadrupoleMoment", value));
    assert(value == "Databases/quadrupoleMoment.txt");
}

void parse_vector() {
    std::vector<std::string> vec;
    std::string section = "  isOn: false\n"
                          "  folder: simulation1\n"
                          "  dataFiles:  \n"
                          "    - eedf\n"
                          "    - swarmParameters  \n"
                          "    -   rateCoefficients\t \n"
                          "    -  \t powerBalance\n"
                          "    - lookUpTable";
    std::string fieldName = "dataFiles";

    std::vector<std::string> correct_vec = {"eedf", "swarmParameters", "rateCoefficients",
                                            "powerBalance", "lookUpTable"};

    assert(Parse::setField(section, fieldName, vec));
    for (uint8_t i = 0; i < vec.size(); ++i) {
        assert(vec[i] == correct_vec[i]);
    }
}

void parse_section() {
    std::string content = "  gasProperties:                       \n"
                          "    mass: Databases/masses.txt\n"
                          "    fraction:\n"
                          "      - N2 = 1\n"
                          "    harmonicFrequency: Databases/harmonicFrequencies.txt\n"
                          "    anharmonicFrequency: Databases/anharmonicFrequencies.txt\n"
                          "    rotationalConstant: Databases/rotationalConstants.txt\n"
                          "    electricQuadrupoleMoment: Databases/quadrupoleMoment.txt\n"
                          "    OPBParameter: Databases/OPBParameter.txt\n"
                          "  stateProperties:                     \n"
                          "    energy:";

    std::string correct_result = "    mass: Databases/masses.txt\n"
                                 "    fraction:\n"
                                 "      - N2 = 1\n"
                                 "    harmonicFrequency: Databases/harmonicFrequencies.txt\n"
                                 "    anharmonicFrequency: Databases/anharmonicFrequencies.txt\n"
                                 "    rotationalConstant: Databases/rotationalConstants.txt\n"
                                 "    electricQuadrupoleMoment: Databases/quadrupoleMoment.txt\n"
                                 "    OPBParameter: Databases/OPBParameter.txt";

    std::string sectionBuffer, sectionTitle = "gasProperties";

    assert(Parse::getSection(content, sectionTitle, sectionBuffer));
    assert(sectionBuffer == correct_result);
}

void parse_tests_phase_one() {
    test(parse_bool);
    test(parse_double);
    test(parse_string);
    test(parse_vector);
    test(parse_section);
}

/** PHASE TWO **/

void parse_value_range() {
    std::string rString = "logspace(11.4574, 23.5, 8)";
    std::string rStringOther = "logspace( 748.,23.5,8)";
    std::string vString = "23.48392";
    double value = 0;

    assert(Parse::getFirstValue(rString, value));
    assert(value == 11.4574);

    assert(Parse::getFirstValue(rStringOther, value));
    assert(value == 748.);

    assert(Parse::getFirstValue(vString, value));
    assert(value == 23.48392);
}

void parse_tests_phase_two() {
    test(parse_value_range);
}
