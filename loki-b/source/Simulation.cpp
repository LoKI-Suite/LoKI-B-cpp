//
// Created by daan on 13-5-19.
//

#include "Simulation.h"
#include <chrono>

namespace loki {
    Simulation::Simulation(const loki::Setup &setup)
            : workingConditions(setup.workingConditions),
              jobManager(&workingConditions) {

        if (setup.electronKinetics.eedfType != EedfType::boltzmann)
        {
            throw std::runtime_error("Only EEDF type 'boltzmann' is supported at present.");
        }

        if (setup.electronKinetics.isOn) {
            initializeJobs(setup.workingConditions);

            electronKinetics = std::make_unique<ElectronKinetics>(setup.electronKinetics, &workingConditions);
            electronKinetics->obtainedNewEedf.addListener(&ResultEvent::emit, &obtainedResults);

            if (setup.output.isOn) {
                output.reset(new Output(setup, &workingConditions, &jobManager));

                electronKinetics->obtainedNewEedf.addListener(&Output::saveCycle, output.get());
                output->simPathExists.addListener(&Event<std::string>::emit, &outputPathExists);
            }
        }
    }
    Simulation::Simulation(const json_type& cnf)
            : workingConditions( cnf.at("workingConditions")),
              jobManager(&workingConditions) {

        if (Enumeration::getEedfType(cnf.at("electronKinetics").at("eedfType")) != EedfType::boltzmann)
        {
            throw std::runtime_error("Only EEDF type 'boltzmann' is supported at present.");
        }

        if (cnf.at("electronKinetics").at("isOn")) {
            initializeJobs(cnf.at("workingConditions"));

            electronKinetics = std::make_unique<ElectronKinetics>(cnf.at("electronKinetics"), &workingConditions);
            electronKinetics->obtainedNewEedf.addListener(&ResultEvent::emit, &obtainedResults);

            if (cnf.at("output").at("isOn")) {
                output.reset(new Output(cnf, &workingConditions, &jobManager));

                electronKinetics->obtainedNewEedf.addListener(&Output::saveCycle, output.get());
                output->simPathExists.addListener(&Event<std::string>::emit, &outputPathExists);
            }
        }
    }

    void Simulation::run() {
        if (electronKinetics.get()) {
            jobManager.prepareFirstJob();
            do {
                electronKinetics->solve();
            } while (jobManager.nextJob());
        }
    }

    Simulation::~Simulation() {
    }

    void Simulation::initializeJobs(const WorkingConditionsSetup &setup) {

        // Repeat this for any other fields that can be declared as a range.
        try {
            const Range range{setup.reducedField};
            jobManager.addJob({"Reduced Field", &WorkingConditions::updateReducedField,range});
        }
        catch(std::exception& exc)
        {
                Log<Message>::Error("Error setting up reduced field: '" + std::string(exc.what()));
        }
    }
    void Simulation::initializeJobs(const json_type &cnf) {

        // Repeat this for any other fields that can be declared as a range.
        try {
            const Range range{cnf.at("reducedField")};
            jobManager.addJob({"Reduced Field", &WorkingConditions::updateReducedField,range});
        }
        catch(std::exception& exc)
        {
                Log<Message>::Error("Error setting up reduced field: '" + std::string(exc.what()));
        }
    }

}
