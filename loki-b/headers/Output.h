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
        const Grid *grid;
        const WorkingConditions *workingConditions;

        std::string folder;

        bool saveEedf{false}, savePower{false}, saveSwarm{false}, saveRates{false}, saveTable{false};

        bool initTable{true};
    public:
        event simPathExists;

    public:
        explicit Output(const OutputSetup &setup, const Grid *grid, const WorkingConditions *workingConditions);

        void saveCycle(const Vector &eedf, const Power &power, const std::vector<EedfGas *> &gasses,
                       const SwarmParameters &swarmParameters,
                       const std::vector<RateCoefficient> &rateCoefficients,
                       const std::vector<RateCoefficient> &extraRateCoefficients, const Vector &firstAnisotropy);

        void plot(const std::string &title, const std::string &xlabel, const std::string &ylabel,
                  const Vector &x, const Vector &y);

    private:
        void writeEedf(const Vector &eedf, const Vector &firstAnisotropy, const Vector &energies);

        void writeSwarm(const SwarmParameters &swarmParameters);

        void writePower(const Power &power, const std::vector<EedfGas *> &gasses);

        void writeRateCoefficients(const std::vector<RateCoefficient> &rateCoefficients,
                                   const std::vector<RateCoefficient> &extraRateCoefficients);

        void writeLookuptable(const Power &power, const SwarmParameters &swarmParameters);
    };
}


#endif //LOKI_CPP_OUTPUT_H
