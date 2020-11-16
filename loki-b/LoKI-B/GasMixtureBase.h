#ifndef LOKI_CPP_GASMIXTUREBASE_H
#define LOKI_CPP_GASMIXTUREBASE_H

#include "LoKI-B/GasBase.h"
#include "LoKI-B/Setup.h"
#include "LoKI-B/WorkingConditions.h"
#include "LoKI-B/json.h"
#include "LoKI-B/Parse.h"
#include "LoKI-B/Log.h"

#include <memory>
#include <vector>
#include <regex>

namespace loki
{

/** This define accepts the name of a gas property (e.g. mass) and will try load the corresponding database
 *  file. If this fails it will output a warning (but nothing more). If the file is succesfully loaded, it
 *  will then try to set the property for all gases in the mixture. This version will not throw an error
 *  in any case, and can therefore be used for gas properties that are not mandatory for every gas (e.g.
 *  atomic gases do not have a harmonicFrequency).
 *
 *  It assumes that the passed property has the same name in the Gas and GasPropertiesSetup structures.
 *  Furthermore, it assumes that a std::string fileBuffer has been declared beforehand.
 *
 *  \todo For now, the electron gas (name=="e") is skipped. Decide how to handle electron properties.
 *  \todo Error reporting should be fixed. The file name is not displayed, the field instead.
 */
template <typename GasListType, typename HandlerType>
void readGasPropertyFile(const GasListType& gasList,
                    const std::string& fileName,
                    const std::string& propertyName,
                    bool required,
                    HandlerType handler)
{
    std::string fileBuffer;
    if (!fileName.empty() && Parse::stringBufferFromFile(fileName, fileBuffer))
    {
        Log<Message>::Notify("Configuring gas property '" + propertyName
              + "', using file '" + fileName + "'.");
        for (auto &gas : gasList)
        {
            if (gas->name == "e")
                continue;
            double value;
            const std::regex r(R"((?:^|\n))" + gas->name + R"(\s+(\S*)\s*)");
            std::smatch m;
            if (std::regex_search(fileBuffer, m, r) && Parse::getValue(m[1],value))
            {
                handler(*gas,value);
            }
            else
            {
                if (required)
                    Log<GasPropertyError>::Error(propertyName + " in gas " + gas->name);
                else
                    Log<GasPropertyError>::Warning(propertyName + " in gas " + gas->name);
            }
        }
    }
    else
    {
        if (required)
            Log<FileError>::Error(propertyName);
        else
            Log<FileError>::Warning(propertyName);
    }
}


/// \todo Error reporting should be fixed. The file name is not displayed, the field instead.
template <typename GasListType, typename HandlerType>
void readGasProperty(const GasListType& gasList,
                    const json_type& cnf,
                    const std::string& propertyName,
                    bool required,
                    HandlerType handler)
{
    const std::string fileName = cnf.contains(propertyName)
        ? cnf.at(propertyName).get<std::string>() : std::string{};
    readGasPropertyFile(gasList,fileName,propertyName,required,handler);
}

/** Gas mixture base class...
 */
class GasMixtureBase
{
  public:
    using Gases = std::vector<std::unique_ptr<GasBase>>;
    virtual ~GasMixtureBase()
    {
    }
    void print(std::ostream &os);
    /** Checks whether the sum of the gas fractions is equal to 1.
     */
    void checkGasFractions();
    void checkPopulations();
    /** Tries to find a gas in the gas mixture specified by a supplied name. If the
     *  gas is present, a pointer to this gas is returned, otherwise a nullptr is
     *  returned.
     */
    GasBase *findGas(const std::string &name);
    /** Tries to find a state in the gas mixture specified by a supplied name. If the
     *  state is present, a pointer to this state is returned, otherwise a nullptr is
     *  returned.
     *
     *  Note that when entry describes a series of states (thus it contains a wildcard
     *  character) then a pointer to the first sibling is returned.
     */
    GasBase::StateBase *findState(const StateEntry &entry);

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
    virtual void loadStateProperties(const StatePropertiesSetup &setup, const WorkingConditions *workingConditions);
    virtual void loadStateProperties(const json_type &cnf, const WorkingConditions *workingConditions);
    /** Loads the gas properties from the database files specified in the main input
     *  file. It is declared virtual such that inherited classes can override it to
     *  declare any extra properties that they might introduce. They can then call
     *  the base class version to load the properties declared in the base class.
     *
     *  \todo It will be easier to implement this in the GasBase class. The
     *        fact that we will then need to open and query the properties file
     *        a few times is no problem, since even in a very complex mixture
     *        there will be only a few gases. The advantage is that we can then
     *        (in the derived mixture class, instead of here) take special measures
     *        for special gases (CAR gases, the electron), which need more or fewer
     *        configuration data. Some properties can then be stored in a property
     *        map that is stored in the (deived) mixture class, instead of in the
     *        GasBase class, where not all properties are always relevant.
     *        This also applies to the other loadGasProperties overload.
     */
    virtual void loadGasProperties(const GasPropertiesSetup &setup);
    virtual void loadGasProperties(const json_type &cnf);

  protected:
    GasBase *addGas(GasBase *gas);

  private:
    Gases m_gases;
};

} // namespace loki

#endif // LOKI_CPP_GASMIXTUREBASE_H
