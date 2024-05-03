#ifndef LOKI_CPP_GASMIXTUREBASE_H
#define LOKI_CPP_GASMIXTUREBASE_H

#include "LoKI-B/Gas.h"
#include "LoKI-B/WorkingConditions.h"
#include "LoKI-B/json.h"
#include "LoKI-B/Parse.h"
#include "LoKI-B/Log.h"

#include <map>
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
 *  \todo Update comments. There is no GasPropertiesSetup structure, only a local fileBuffer.
 */
template <typename GasListType, typename HandlerType>
void readGasPropertyFile(const std::filesystem::path &basePath, 
                    const GasListType& gasList,
                    const std::string& fileName,
                    const std::string& propertyName,
                    bool required,
                    HandlerType handler)
{
    std::filesystem::path path(fileName);

    if (path.is_relative()) {
        path = basePath.parent_path() / path;
    }

    std::string fileBuffer;
    if (!path.empty() && Parse::stringBufferFromFile(path, fileBuffer))
    {
        Log<Message>::Notify("Configuring gas property '", propertyName,
              "', using file '", path,  "'.");
        for (auto &gas : gasList)
        {
            if (gas->name() == "e")
                continue;
            double value;
            const std::regex r(R"((?:^|\n))" + gas->name() + R"(\s+(\S*)\s*)");
            std::smatch m;
            if (std::regex_search(fileBuffer, m, r) && Parse::getValue(m[1],value))
            {
                handler(*gas,value);
            }
            else
            {
                if (required)
                    Log<GasPropertyError>::Error(propertyName + " in gas " + gas->name());
                else
                    Log<GasPropertyError>::Warning(propertyName + " in gas " + gas->name());
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

/** This define accepts the name of a gas property (e.g. mass) and will try to parse the json object.
 *  If this fails it will output a warning (but nothing more). If the json is succesfully loaded, it
 *  will then try to set the property for all gases in the mixture. This version will not throw an error
 *  in any case, and can therefore be used for gas properties that are not mandatory for every gas (e.g.
 *  atomic gases do not have a harmonicFrequency).
 *
 *  It assumes that the passed property has the same name in the Gas and GasPropertiesSetup structures.
 *  Furthermore, it assumes that a std::string fileBuffer has been declared beforehand.
 *
 *  \todo For now, the electron gas (name=="e") is skipped. Decide how to handle electron properties.
 */
template <typename GasListType, typename HandlerType>
void readGasPropertyJson(const GasListType& gasList,
                    const json_type& cnf,
                    const std::string& propertyName,
                    bool required,
                    HandlerType handler)
{

    Log<Message>::Notify("Configuring gas property '" + propertyName + "'.");
    for (auto &gas : gasList)
    {
        if (gas->name() == "e")
            continue;
	bool found = false;
        if (cnf.contains(gas->name()))
	{
            found = true;
            handler(*gas, cnf.at(gas->name()).template get<double>());
        }
        if (!found)
        {
            if (required)
                Log<GasPropertyError>::Error(propertyName + " in gas " + gas->name());
            else
                Log<GasPropertyError>::Warning(propertyName + " in gas " + gas->name());
        }
    }
}

/// \todo Error reporting should be fixed. The file name is not displayed, the field instead.
template <typename GasListType, typename HandlerType>
void readGasProperty(const std::filesystem::path &basePath,
                     const GasListType& gasList,
                     const json_type& cnf,
                     const std::string& propertyName,
                     bool required,
                     HandlerType handler)
{
    if (cnf.contains(propertyName) && cnf.at(propertyName).is_object())
    {
        readGasPropertyJson(gasList,cnf.at(propertyName),propertyName,required,handler);
    }
    else
    {
        const auto fileName = cnf.contains(propertyName)
            ? cnf.at(propertyName).get<std::string>() : std::string{};
        readGasPropertyFile(basePath,gasList,fileName,propertyName,required,handler);
    }
}

/** Gas mixture base class...
 */
class GasMixture
{
  public:
    using Gases = std::vector<std::unique_ptr<Gas>>;
    ~GasMixture();
    void print(std::ostream &os);
    /** Checks whether the sum of the gas fractions is equal to 1.
     */
    void checkGasFractions();
    void checkPopulations();
    // Vector of pointers to all the gases in the mixture.
    const Gases &gases() const { return m_gases; }
    /** Tries to add a new gas based on a given name and returns a pointer to it. If a
     *  gas with the same name already exists it returns a pointer to that gas.
     */
    Gas *ensureGas(const std::string &name);
    /** Tries to find a gas in the gas mixture specified by a supplied name. If the
     *  gas is present, a pointer to this gas is returned, otherwise a nullptr is
     *  returned.
     */
    Gas *findGas(const std::string &name);
    /** Tries to find a state in the gas mixture specified by a supplied name. If the
     *  state is present, a pointer to this state is returned, otherwise a nullptr is
     *  returned.
     *
     *  Note that when entry describes a series of states (thus it contains a wildcard
     *  character) then a pointer to the first sibling is returned.
     */
    Gas::State *findState(const StateEntry &entry);
    /** Returns a container of matches. Empty (no match), one state (no owildcard)
     *  or a vector with one or more state pointers (entry has a wildcard).
     *  \todo reconsider all the state accessors. Which interfaces do we really need.
     */
    Gas::State::ChildContainer findStates(const StateEntry &entry);
    /** Loads the data concerning a single property of the states (energy, statistical
     *  weight, or population) from a vector of entries as supplied by the setup object.
     *  First it determines whether the current entry requires loading by direct value,
     *  file or function, then it acts accordingly. The property that needs to be set
     *  is defined by the propertyType variable. Furthermore, it also needs a pointer
     *  to the WorkingConditions structure to access the argument map (which maps its
     *  member variables to names by which they are addressed in the input files).
     */
    void loadStateProperty(const std::filesystem::path &basePath, const std::vector<std::string> &entryVector, 
                           StatePropertyType propertyType, const WorkingConditions *workingConditions);
    /** Evaluates the densities of all states in the state tree. \todo Explain. All states of all gases?
     */
    void evaluateReducedDensities();
    /** Loads the state properties as specified in the input file. It also calls
     *  checkPopulations.
     */
    void loadStateProperties(const std::filesystem::path &basePath, const json_type &cnf, 
                             const WorkingConditions *workingConditions);
    /** Loads the gas properties from the database files specified in the main input
     *  file.
     *
     *  \todo It will be easier to implement this in the Gas class. The
     *        fact that we will then need to open and query the properties file
     *        a few times is no problem, since even in a very complex mixture
     *        there will be only a few gases. The advantage is that we can then
     *        (in the derived mixture class, instead of here) take special measures
     *        for special gases (CAR gases, the electron), which need more or fewer
     *        configuration data. Some properties can then be stored in a property
     *        map that is stored in the (deived) mixture class, instead of in the
     *        Gas class, where not all properties are always relevant.
     */
    void loadGasProperties(const std::filesystem::path &basePath, const json_type &cnf);

    using StateMap = std::map<std::string, Gas::State *>;
    /** This overload accepts only a reference to a StateEntry object. This is the main
     *  function to call when a new state is to be added to the mixture. It will add the
     *  corresponding gas and ancestor states if the do not yet exist, before adding the
     *  state itself. If the state already exists, then a pointer to this state is
     *  returned instead.
     *
     *  Note that this container only contains states that are actually present in the
     *  mixture, not their parents that are crated implicitly.
     */
    Gas::State *ensureState(const StateEntry &entry);
    Gas::State *findStateById(const std::string &stateId);
  private:
    Gas *addGas(const std::string& name);
    Gases m_gases;
    StateMap m_states;
};

} // namespace loki

#endif // LOKI_CPP_GASMIXTUREBASE_H
