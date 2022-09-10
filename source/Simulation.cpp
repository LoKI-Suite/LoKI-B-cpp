/** \file
 *
 *  Implementation of LoKI-B's Simulation class.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2022 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
 *  M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
 *  L. L. Alves
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  \author Daan Boer and Jan van Dijk (C++ version)
 *  \date   13 May 2019 (first C++ version)
 */

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

        m_electronKinetics = std::make_unique<ElectronKineticsBoltzmann>(setup.electronKinetics, &m_workingConditions);
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

        m_electronKinetics = std::make_unique<ElectronKineticsBoltzmann>(cnf.at("electronKinetics"), &m_workingConditions);
        m_electronKinetics->obtainedNewEedf.addListener(&ResultEvent::emit, &m_obtainedResults);
    }
    Log<Message>::Notify("Simulation has been set up", ", number of parameters = ", m_jobManager.dimension(),
                         ", number of jobs = ", m_jobManager.njobs());
}

void Simulation::run()
{
    if (m_electronKinetics.get())
    {
        m_jobManager.prepareFirstJob();
        do
        {
            m_workingConditions.setCurrentJobFolder(m_jobManager.getCurrentJobFolder());
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
            "ReducedField",
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
            "ReducedField",
            std::bind(&WorkingConditions::updateReducedField, std::ref(m_workingConditions), std::placeholders::_1),
            Range::create(cnf.at("reducedField")));
    }
    catch (std::exception &exc)
    {
        Log<Message>::Error("Error setting up reduced field: '" + std::string(exc.what()));
    }
}

} // namespace loki
