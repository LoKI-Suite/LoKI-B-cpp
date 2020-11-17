//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_WORKINGCONDITIONS_H
#define LOKI_CPP_WORKINGCONDITIONS_H

#include "LoKI-B/Event.h"
#include "LoKI-B/Setup.h"
#include "LoKI-B/json.h"
#include "LoKI-B/Constant.h"
#include <string>
#include <map>

namespace loki
{

/** The WorkingConditions class stores all global variables
 *  that make up the working environment of the simulation.
 *  It is also responsible for propagating any change in the
 *  working conditions to the relevant classes.
 */
class WorkingConditions
{
  public:

    // constructors and destructor

    WorkingConditions(const WorkingConditionsSetup &setup);
    WorkingConditions(const json_type &cnf);
    ~WorkingConditions() = default;
    // Copying this object is not allowed.
    WorkingConditions(const WorkingConditions &other) = delete;

    double gasTemperature() const { return m_gasTemperature; }
    double gasDensity() const { return m_gasDensity; }
    double electronDensity() const { return m_electronDensity; }
    double reducedField() const { return m_reducedField; }
    double reducedFieldSI() const { return m_reducedField * Constant::SI::Townsend; }
    double excitationFrequency() const { return m_excitationFrequency; }
    double reducedExcFreqSI() const { return m_excitationFrequency * 2 * Constant::pi / m_gasDensity; }

    using ArgumentMap = std::map<std::string, const double *>;
    const ArgumentMap& argumentMap() const { return m_argumentMap; }

    void updateReducedField(double value);
    void updateElectronTemperature(double value);

    /** \todo Technically we only need one 'updatedGasTemperature' event,
     *        since we can control the order of the listeners.
     */
    /** \todo Document the individual events carefully. What are the intent
     *        and expectations?
     */
    /** \todo Check which events are actually used. Note that this class has
     *        only two 'update' members, for the reduced field and for the
     *        electron temperature. Update: only those two events are used
     *        at present.
     */
    Event<> updatedReducedField;
    Event<> updatedElectronTemperature;
    //Event<> updatedGasPressure;
    //Event<> updatedGasTemperature1;
    //Event<> updatedGasTemperature2;
    //Event<> updatedGasDensity;
    //Event<> updatedElectronDensity;
    //Event<> updatedChamberLength;
    //Event<> updatedExcitationFrequency;
private:
    const double m_gasPressure;
    const double m_gasTemperature;
    const double m_gasDensity;
    double m_electronDensity;
    /** \todo It appears that electronTemperature is not used anywhere
     *  at present, maybe because only EEDF type 'boltzmann' has been
     *  implemented? It is not used, and set only by
     *  evaluateSwarmParameters() at present. Check this.
     */
    double m_electronTemperature;
    /// \todo chamberLength is not used anywhere at present
    //double chamberLength;
    /// \todo chamberRadius is not used anywhere at present
    //double chamberRadius;
    double m_reducedField;
    double m_excitationFrequency;
    ArgumentMap m_argumentMap;
};

} // namespace loki

#endif // LOKI_CPP_WORKINGCONDITIONS_H
