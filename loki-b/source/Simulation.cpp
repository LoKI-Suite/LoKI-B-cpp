//
// Created by daan on 13-5-19.
//

#include "Simulation.h"
#include <chrono>

namespace loki {
    Simulation::Simulation(const loki::Setup &setup)
            : workingConditions(setup.workingConditions, setup.electronKinetics.eedfType),
              enableKinetics(setup.electronKinetics.isOn), enableOutput(setup.output.isOn),
              jobManager(&workingConditions) {

        if (enableKinetics) {
            initializeJobs(setup.workingConditions);

            electronKinetics = std::make_unique<ElectronKinetics>(setup.electronKinetics, &workingConditions);
            electronKinetics->obtainedNewEedf.addListener(&ResultEvent::emit, &obtainedResults);

            if (enableOutput) {
                output = new Output(setup, &workingConditions, &jobManager);

                electronKinetics->obtainedNewEedf.addListener(&Output::saveCycle, output);
                output->simPathExists.addListener(&Event<std::string>::emit, &outputPathExists);
            }
        }
    }
    Simulation::Simulation(const json_type& cnf)
            : workingConditions(
                  cnf.at("workingConditions"),
                  Enumeration::getEedfType(cnf.at("electronKinetics").at("eedfType")) ),
              enableKinetics(cnf.at("electronKinetics").at("isOn")),
              enableOutput(cnf.at("output").at("isOn")),
              jobManager(&workingConditions) {

        if (enableKinetics) {
            initializeJobs(cnf.at("workingConditions"));

            electronKinetics = std::make_unique<ElectronKinetics>(cnf.at("electronKinetics"), &workingConditions);
            electronKinetics->obtainedNewEedf.addListener(&ResultEvent::emit, &obtainedResults);

            if (enableOutput) {
                output = new Output(cnf, &workingConditions, &jobManager);

                electronKinetics->obtainedNewEedf.addListener(&Output::saveCycle, output);
                output->simPathExists.addListener(&Event<std::string>::emit, &outputPathExists);
            }
        }
    }

    void Simulation::run() {
        if (enableKinetics) {
            if (multipleSimulations) {
                do {
                    electronKinetics->solve();
                } while (jobManager.nextJob());
            } else {
                electronKinetics->solve();
            }

        }
    }

    Simulation::~Simulation() {
        if (enableOutput) delete output;
    }

    void Simulation::initializeJobs(const WorkingConditionsSetup &setup) {

        // Repeat this statement for any other fields that can be declared as a range.
        if (!Parse::isNumerical(setup.reducedField)) {
            if (!initializeJob("Reduced Field", setup.reducedField, &WorkingConditions::updateReducedField)) {
                Log<Message>::Error("Reduced field entry in input file is ill formatted.");
            }
        }
    }
    void Simulation::initializeJobs(const json_type &cnf) {

        // Repeat this statement for any other fields that can be declared as a range.
        if (cnf.at("reducedField").type()==json_type::value_t::string) {
            if (!initializeJob("Reduced Field", cnf.at("reducedField"), &WorkingConditions::updateReducedField)) {
                Log<Message>::Error("Reduced field entry in input file is ill formatted.");
            }
        }
    }

    bool Simulation::initializeJob(const std::string &name, const std::string &valueString,
                                   void (WorkingConditions::*callback)(double)) {
        bool success = false;

        Range range = Parse::getRange(valueString, success);

        if (!success) return false;

        jobManager.addJob({name, &WorkingConditions::updateReducedField, range});
        multipleSimulations = true;

        return true;
    }
}
