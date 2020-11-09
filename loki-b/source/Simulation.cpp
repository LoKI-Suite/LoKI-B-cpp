//
// Created by daan on 13-5-19.
//

#include "LoKI-B/Simulation.h"

namespace loki
{

Simulation::Simulation(const loki::Setup &setup)
    : m_workingConditions(setup.workingConditions),
    m_jobManager()
{
    if (setup.electronKinetics.eedfType != EedfType::boltzmann)
    {
        throw std::runtime_error("Only EEDF type 'boltzmann' is supported at present.");
    }

    if (setup.electronKinetics.isOn)
    {
        initializeJobs(setup.workingConditions);

        m_electronKinetics = std::make_unique<ElectronKinetics>(setup.electronKinetics, &m_workingConditions);
        m_electronKinetics->obtainedNewEedf.addListener(&ResultEvent::emit, &m_obtainedResults);
    }
    Log<Message>::Notify("Simulation has been set up", ", number of parameters = ", m_jobManager.dimension(),
                         ", number of jobs = ", m_jobManager.njobs());
}

Simulation::Simulation(const json_type &cnf)
    : m_workingConditions(cnf.at("workingConditions")),
    m_jobManager()
{
    if (getEedfType(cnf.at("electronKinetics").at("eedfType")) != EedfType::boltzmann)
    {
        throw std::runtime_error("Only EEDF type 'boltzmann' is supported at present.");
    }

    if (cnf.at("electronKinetics").at("isOn"))
    {
        initializeJobs(cnf.at("workingConditions"));

        m_electronKinetics = std::make_unique<ElectronKinetics>(cnf.at("electronKinetics"), &m_workingConditions);
        m_electronKinetics->obtainedNewEedf.addListener(&ResultEvent::emit, &m_obtainedResults);
    }
    Log<Message>::Notify("Simulation has been set up", ", number of parameters = ", m_jobManager.dimension(),
                         ", number of jobs = ", m_jobManager.njobs());
}

void Simulation::configureOutput(Output* output)
{
    m_output.reset(output);
    if (m_electronKinetics.get())
    {
        m_electronKinetics->obtainedNewEedf.addListener(&Output::saveCycle, m_output.get());
    }
}

void Simulation::run()
{
    if (m_electronKinetics.get())
    {
        m_jobManager.prepareFirstJob();
        do
        {
            m_electronKinetics->solve();
        } while (m_jobManager.prepareNextJob());
    }
}

Simulation::~Simulation()
{
}

void Simulation::initializeJobs(const WorkingConditionsSetup &setup)
{

    // Repeat this for any other fields that can be declared as a range.
    try
    {
        m_jobManager.addParameter(
            "Reduced Field",
            std::bind(&WorkingConditions::updateReducedField, std::ref(m_workingConditions), std::placeholders::_1),
            Range::create(setup.reducedField));
    }
    catch (std::exception &exc)
    {
        Log<Message>::Error("Error setting up reduced field: '" + std::string(exc.what()));
    }
}

void Simulation::initializeJobs(const json_type &cnf)
{

    // Repeat this for any other fields that can be declared as a range.
    try
    {
        m_jobManager.addParameter(
            "Reduced Field",
            std::bind(&WorkingConditions::updateReducedField, std::ref(m_workingConditions), std::placeholders::_1),
            Range::create(cnf.at("reducedField")));
    }
    catch (std::exception &exc)
    {
        Log<Message>::Error("Error setting up reduced field: '" + std::string(exc.what()));
    }
}

} // namespace loki
