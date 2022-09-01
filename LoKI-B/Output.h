//
// Created by daan on 04-07-2019.
//

#ifndef LOKI_CPP_OUTPUT_H
#define LOKI_CPP_OUTPUT_H

#include <string>
#include <functional>

#include "LoKI-B/Grid.h"
#include "LoKI-B/JobSystem.h"
#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/MacroscopicQuantities.h"
#include "LoKI-B/Power.h"
#include "LoKI-B/Setup.h"
#include "LoKI-B/WorkingConditions.h"
#include "LoKI-B/json.h"
#include "LoKI-B/Exports.h"

namespace loki
{

class EedfCollisionDataMixture;

class lokib_export Output
{
public:
    virtual ~Output();
    void saveCycle(const Grid &energyGrid, const Vector &eedf, const WorkingConditions &wc, const Power &power,
                   const EedfCollisionDataMixture& collData, const SwarmParameters &swarmParameters,
                   const Vector &firstAnisotropy);
protected:
    Output(const Setup &setup, const WorkingConditions *workingConditions);
    Output(const json_type &cnf, const WorkingConditions *workingConditions);
    virtual void setDestination(const std::string& subFolder)=0;
    virtual void writeEedf(const Vector &eedf, const Vector &firstAnisotropy, const Vector &energies) const=0;
    virtual void writeSwarm(const SwarmParameters &swarmParameters) const=0;
    virtual void writePower(const Power &power, const EedfCollisionDataMixture &collData) const=0;
    virtual void writeRateCoefficients(const std::vector<RateCoefficient> &rateCoefficients,
                               const std::vector<RateCoefficient> &extraRateCoefficients) const=0;
    virtual void writeLookuptable(const Power &power, const SwarmParameters &swarmParameters) const=0;
    const WorkingConditions *workingConditions;
private:
    const JobManager *jobManager;
    bool saveEedf, savePower, saveSwarm, saveRates, saveTable;
};

class lokib_export FileOutput : public Output
{
public:
    using PathExistsHandler = std::function<void(std::string&)>;
    FileOutput(const Setup &setup, const WorkingConditions *workingConditions,
        const PathExistsHandler& handler);
    FileOutput(const json_type &cnf, const WorkingConditions *workingConditions,
        const PathExistsHandler& handler);
protected:
    virtual void setDestination(const std::string& subFolder);
    virtual void writeEedf(const Vector &eedf, const Vector &firstAnisotropy, const Vector &energies) const;
    virtual void writeSwarm(const SwarmParameters &swarmParameters) const;
    virtual void writePower(const Power &power, const EedfCollisionDataMixture &collData) const;
    virtual void writeRateCoefficients(const std::vector<RateCoefficient> &rateCoefficients,
                               const std::vector<RateCoefficient> &extraRateCoefficients) const;
    virtual void writeLookuptable(const Power &power, const SwarmParameters &swarmParameters) const;
private:
    void createPath(const PathExistsHandler& handler);
    void writeTerm(std::ostream& os, const std::string& name, const std::string& unit, double value, bool plus=false) const;
    std::string m_folder;
    std::string m_subFolder;
    mutable bool m_initTable;
};

class lokib_export JsonOutput : public Output
{
public:
    JsonOutput(json_type& root, const Setup &setup, const WorkingConditions *workingConditions);
    JsonOutput(json_type& root, const json_type &cnf, const WorkingConditions *workingConditions);
protected:
    virtual void setDestination(const std::string& subFolder);
    virtual void writeEedf(const Vector &eedf, const Vector &firstAnisotropy, const Vector &energies) const;
    virtual void writeSwarm(const SwarmParameters &swarmParameters) const;
    virtual void writePower(const Power &power, const EedfCollisionDataMixture& collData) const;
    virtual void writeRateCoefficients(const std::vector<RateCoefficient> &rateCoefficients,
                               const std::vector<RateCoefficient> &extraRateCoefficients) const;
    virtual void writeLookuptable(const Power &power, const SwarmParameters &swarmParameters) const;
private:
    /** This produces a member of the form: name: { "value": value, "unit": unit }.
     *  As an example, makeQuantity("Te", 2.0, "eV") produces and returns an object that contains
     *  \verbatim
          "Te": { "value": = 2.0, "unit": "eV" }. \endverbatim
     *
     */
    static json_type makeQuantity(const std::string& name, double value, const std::string unit);
    json_type& m_root;
    json_type* m_active;
};

} // namespace loki

#endif // LOKI_CPP_OUTPUT_H
