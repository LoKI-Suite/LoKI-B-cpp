#ifndef LOKI_CPP_GASMIXTUREBASE_H
#define LOKI_CPP_GASMIXTUREBASE_H

#include "GasBase.h"
#include "Setup.h"
#include "json.h"
#include "WorkingConditions.h"

#include <vector>
#include <memory>

/* -- GAS_PROPERTY --
 * This define accepts the name of a gas property (e.g. mass) and will try load the corresponding database
 * file. If this fails it will output a warning (but nothing more). If the file is succesfully loaded, it
 * will then try to set the property for all gases in the mixture. This version will not throw an error
 * in any case, and can therefore be used for gas properties that are not mandatory for every gas (e.g.
 * atomic gases do not have a harmonicFrequency).
 *
 * It assumes that the passed property has the same name in the Gas and GasPropertiesSetup structures.
 * Furthermore, it assumes that a std::string fileBuffer has been declared beforehand.
 */

#define GAS_PROPERTY_TEMPL(GASLIST,PROPERTY,SEVERITY) \
{ \
    std::string fileBuffer; \
    std::cout << "Configuring gas property '" << #PROPERTY << "', using file '" << setup.PROPERTY << "'." << std::endl; \
    if (Parse::stringBufferFromFile(setup.PROPERTY, fileBuffer)) { \
        for (auto& gas : GASLIST) { \
            if(!Parse::gasProperty(gas->name, gas->PROPERTY, fileBuffer)) { \
                Log<GasPropertyError>::SEVERITY(#PROPERTY " in gas " + gas->name); \
            } \
        } \
    } \
    else { \
        Log<FileError>::SEVERITY(setup.PROPERTY); \
    } \
}

#define GAS_PROPERTY(GASLIST,PROPERTY) GAS_PROPERTY_TEMPL(GASLIST,PROPERTY,Warning)
#define R_GAS_PROPERTY(GASLIST,PROPERTY) GAS_PROPERTY_TEMPL(GASLIST,PROPERTY,Error)

// JSON versions. The filename is obtained from cnf, the JSON version of the input file
/// \todo Error reporting should be fixed. The file name is not displayed, the field instead.
#define GAS_PROPERTY_JSON_TEMPL(GASLIST,cnf_object,PROPERTY,SEVERITY) \
{ \
    std::string fileBuffer; \
    if (cnf_object.contains(#PROPERTY) && Parse::stringBufferFromFile(cnf_object.at(#PROPERTY), fileBuffer)) { \
        std::cout << "Configuring gas property '" << #PROPERTY << "', using file '" << cnf_object.at(#PROPERTY) << "'." << std::endl; \
        for (auto& gas : GASLIST) { \
            if(!Parse::gasProperty(gas->name, gas->PROPERTY, fileBuffer)) { \
                Log<GasPropertyError>::SEVERITY(#PROPERTY " in gas " + gas->name); \
            } \
        } \
    } \
    else { \
        Log<FileError>::SEVERITY(#PROPERTY); \
    } \
}

#define GAS_PROPERTY_JSON(GASLIST,cnf_object,PROPERTY) GAS_PROPERTY_JSON_TEMPL(GASLIST,cnf_object,PROPERTY,Warning);
#define R_GAS_PROPERTY_JSON(GASLIST,cnf_object,PROPERTY) GAS_PROPERTY_JSON_TEMPL(GASLIST,cnf_object,PROPERTY,Error);

namespace loki {

    /** Gas mixture base class...
     */
    class GasMixtureBase
    {
    public:
        using Gases = std::vector<std::unique_ptr<GasBase>>;
        virtual ~GasMixtureBase(){}
        GasBase* addGas(GasBase* gas);
        void print(std::ostream& os);
        /** Checks whether the sum of the gas fractions is equal to 1.
         */
        void checkGasFractions();
        void checkPopulations();
        /** Tries to find a gas in the gas mixture specified by a supplied name. If the
         *  gas is present, a pointer to this gas is returned, otherwise a nullptr is
         *  returned.
         */
        GasBase* findGas(const std::string &name);
        /** Tries to find a state in the gas mixture specified by a supplied name. If the
         *  state is present, a pointer to this state is returned, otherwise a nullptr is
         *  returned.
         *
         *  Note that when entry describes a series of states (thus it contains a wildcard
         *  character) then a pointer to the first sibling is returned.
         */
        GasBase::StateBase* findState(const StateEntry &entry);

        /** Loads the data concerning a single property of the states (energy, statistical
         *  weight, or population) from a vector of entries as supplied by the setup object.
         *  First it determines whether the current entry requires loading by direct value,
         *  file or function, then it acts accordingly. The property that needs to be set
         *  is defined by the propertyType variable. Furthermore, it also needs a pointer
         *  to the WorkingConditions structure to access the argument map (which maps its
         *  member variables to names by which they are addressed in the input files).
         */
        void loadStateProperty(const std::vector<std::string> &entryVector, StatePropertyType propertyType,
                               const WorkingConditions *workingConditions);
        /** Evaluates the densities of all states in the state tree. \todo Explain. All states of all gases?
         */
        void evaluateStateDensities();
        /** Loads the state properties as specified in the input file. It also calls
         *  checkPopulations.
         */
        virtual void loadStateProperties(const StatePropertiesSetup &setup,
                                         const WorkingConditions *workingConditions);
        virtual void loadStateProperties(const json_type &cnf,
                                         const WorkingConditions *workingConditions);
        /** Loads the gas properties from the database files specified in the main input
         *  file. It is declared virtual such that inherited classes can override it to
         *  declare any extra properties that they might introduce. They can then call
         *  the base class version to load the properties declared in the base class.
         */
        virtual void loadGasProperties(const GasPropertiesSetup &setup);
        virtual void loadGasProperties(const json_type &cnf);
  private:
        Gases m_gases;
    };



} // namespace loki

#endif // LOKI_CPP_GASMIXTUREBASE_H

