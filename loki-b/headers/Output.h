//
// Created by daan on 04-07-2019.
//

#ifndef LOKI_CPP_OUTPUT_H
#define LOKI_CPP_OUTPUT_H

#include <string>

#include "Setup.h"
#include "MacroscopicQuantities.h"
#include "Power.h"
#include "LinearAlgebra.h"
#include "Grid.h"
#include "Event.h"
#include "WorkingConditions.h"
#include "EedfGas.h"

namespace loki {
    class Output {
        const WorkingConditions *workingConditions;

        std::string folder, subFolder;

        const JobManager *jobManager;

        bool saveEedf{false}, savePower{false}, saveSwarm{false}, saveRates{false}, saveTable{false};

        bool initTable{true};

        std::string inputFile;
    public:
        Event <std::string> simPathExists;

    public:
        explicit Output(const Setup &setup, const WorkingConditions *workingConditions, const JobManager *jobManager);

        void createPath();

        void saveCycle(const Grid &energyGrid, const Vector &eedf, const WorkingConditions &wc, const Power &power,
                       const std::vector<EedfGas *> &gasses, const SwarmParameters &swarmParameters,
                       const std::vector<RateCoefficient> &rateCoefficients,
                       const std::vector<RateCoefficient> &extraRateCoefficients, const Vector &firstAnisotropy);

    private:
        void writeInputFile();

        void writeEedf(const Vector &eedf, const Vector &firstAnisotropy, const Vector &energies);

        void writeSwarm(const SwarmParameters &swarmParameters);

        void writePower(const Power &power, const std::vector<EedfGas *> &gasses);

        void writeRateCoefficients(const std::vector<RateCoefficient> &rateCoefficients,
                                   const std::vector<RateCoefficient> &extraRateCoefficients);

        void writeLookuptable(const Power &power, const SwarmParameters &swarmParameters);
    };
}


#endif //LOKI_CPP_OUTPUT_H
