//
// Created by daan on 2-5-19.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>

#include "Setup.h"

/*
 * The SET definition provides a shorthand for setting a member variable of a setup struct.
 * The user has to pass
 *
 * The '#property' yields the passed variable name as a string. Hence, SET requires that the
 * designated member variable has the same name as the corresponding field in the input file.
 */
#define SET(section, property) Parse::setField(section, #property, property)

/*
 * The R_SET variant will cause a boolean function to return false when the setField function
 * is unsuccessful. Use this variant for fields that are required to perform the simulation.
 */
#define R_SET(section, property) if (!Parse::setField(section, #property, property)) return false

/*
 * The SUB_STRUCT definition provides a shorthand for setting a sub structure of a SetupBase
 * object.
 */
#define SUB_STRUCT(section, subStruct) parseSubStructure(section, #subStruct, subStruct)

/*
 * The R_SUB_STRUCT requires the sub structures to be present in the input file, and parsed
 * successfully. Use this variant for sub structures that are obligatory.
 */
#define R_SUB_STRUCT(section, subStruct) if (!parseSubStructure(section, #subStruct, subStruct)) return false

namespace loki {

    /*
     * The parseSubStructure function fills a substructure with the data available in the
     * supplied (section of the) input file. This substructure is derived from the
     * SetupBase struct and it is passed by reference. Note that every such derived class
     * has overridden the 'parse' function to fill its fields and substructures as
     * desired.
     */

    bool SetupBase::parseSubStructure(const std::string &content, const std::string &fieldName,
                                      SetupBase &subStruct) {
        std::string sectionContent;

        if (Parse::getSection(content, fieldName, sectionContent)) {
            if (!subStruct.parse(sectionContent)) {
                std::cerr << "[warning] Could not properly parse the " << fieldName
                          << " section. Please check for missing values and improper "
                             "indentation"
                          << std::endl;
                return false;
            }
        } else {
            std::cerr << "[warning] The input file does not contain the "
                      << fieldName << " section." << std::endl;
            return false;
        }

        return true;
    }

    /*
     * The parseFile function is only available to the main Setup class. The user
     * passes the name of the input file. The function will then extract the text
     * from the input file, remove any comments and pass it to the parse function.
     */

    bool Setup::parseFile(const std::string &fileName) {
        std::ifstream file(inputPath + '/' + fileName);

        if (!file.is_open()) {
            std::cerr << "[error] Could not find/open specified file" << std::endl;
            return false;
        }

        // Load the contents of the file into a string buffer.
        std::stringstream stringBuffer;
        stringBuffer << file.rdbuf();

        // Store the file contents in a string and remove any comments.
        std::string fileContent = Parse::removeComments(stringBuffer.str());

        return this->parse(fileContent);
    }

    /*
     * Every derived class of the BaseSetup struct needs to override the parse function.
     * Inside this function the user should call Parse::setField on all the class member
     * variables, and SetupBase::parseSubStructure on all its substructures.
     *
     * At the top of this file, there are some defines that ease this process. However,
     * they do require the member variable to have the same name as specified in the
     * input file.
     *
     * Note that due to the lack of reflection/inspection in C++, we cannot reference
     * members by a string representation or loop through them. Therefore, we have to
     * call either (R_)SET() or (R_)SUB_STRUCT() for each member individually.
     */

    bool Setup::parse(const std::string &sectionContent) {

        R_SUB_STRUCT(sectionContent, workingConditions);
        R_SUB_STRUCT(sectionContent, electronKinetics);
        R_SUB_STRUCT(sectionContent, output);

        return true;
    }

    bool WorkingConditionsSetup::parse(const std::string &sectionContent) {
        // TODO: Check whether 'reducedField' is present in the case that electronKinetics
        //  is enabled (and subsequently that 'electronTemperature' is present when it is
        //  disabled.

        SET(sectionContent, reducedField);
        SET(sectionContent, electronTemperature);
        R_SET(sectionContent, excitationFrequency);
        R_SET(sectionContent, gasPressure);
        R_SET(sectionContent, gasTemperature);
        R_SET(sectionContent, electronDensity);
        R_SET(sectionContent, chamberLength);
        R_SET(sectionContent, chamberRadius);

        return true;
    }

    bool ElectronKineticsSetup::parse(const std::string &sectionContent) {
        R_SET(sectionContent, isOn);
        R_SET(sectionContent, eedfType);
        SET(sectionContent, shapeParameter);
        R_SET(sectionContent, ionizationOperatorType);
        R_SET(sectionContent, growthModelType);
        R_SET(sectionContent, includeEECollisions);
        R_SET(sectionContent, LXCatFiles);
        SET(sectionContent, LXCatFilesExtra);
        SET(sectionContent, effectiveCrossSectionPopulations);
        SET(sectionContent, CARgases);

        R_SUB_STRUCT(sectionContent, gasProperties);
        R_SUB_STRUCT(sectionContent, stateProperties);
        R_SUB_STRUCT(sectionContent, numerics);

        return true;
    }

    bool GasPropertiesSetup::parse(const std::string &sectionContent) {
        R_SET(sectionContent, mass);
        R_SET(sectionContent, fraction);
        R_SET(sectionContent, harmonicFrequency);
        R_SET(sectionContent, anharmonicFrequency);
        R_SET(sectionContent, rotationalConstant);
        R_SET(sectionContent, electricQuadrupoleMoment);
        R_SET(sectionContent, OPBParameter);

        return true;
    }

    bool StatePropertiesSetup::parse(const std::string &sectionContent) {
        R_SET(sectionContent, energy);
        R_SET(sectionContent, statisticalWeight);
        R_SET(sectionContent, population);

        return true;
    }

    bool NumericsSetup::parse(const std::string &sectionContent) {
        SET(sectionContent, maxPowerBalanceRelError);

        R_SUB_STRUCT(sectionContent, energyGrid);
        R_SUB_STRUCT(sectionContent, nonLinearRoutines);

        return true;
    }

    bool EnergyGridSetup::parse(const std::string &sectionContent) {
        R_SET(sectionContent, maxEnergy);
        R_SET(sectionContent, cellNumber);

        SUB_STRUCT(sectionContent, smartGrid);

        return true;
    }

    bool SmartGridSetup::parse(const std::string &sectionContent) {
        R_SET(sectionContent, minEedfDecay);
        R_SET(sectionContent, maxEedfDecay);
        R_SET(sectionContent, updateFactor);

        return true;
    }

    bool NonLinearRoutinesSetup::parse(const std::string &sectionContent) {
        R_SET(sectionContent, algorithm);
        R_SET(sectionContent, mixingParameter);
        R_SET(sectionContent, maxEedfRelError);

        SUB_STRUCT(sectionContent, odeSetParameters);

        return true;
    }

    bool OdeSetParametersSetup::parse(const std::string &sectionContent) {
        R_SET(sectionContent, maxStep);

        return true;
    }

    bool OutputSetup::parse(const std::string &sectionContent) {
        R_SET(sectionContent, isOn);
        R_SET(sectionContent, folder);
        R_SET(sectionContent, dataFiles);

        return true;
    }
}