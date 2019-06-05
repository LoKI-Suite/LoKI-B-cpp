//
// Created by daan on 13-5-19.
//

#include "ElectronKinetics.h"
#include <chrono>

namespace loki {
    ElectronKinetics::ElectronKinetics(const ElectronKineticsSetup &setup, const WorkingConditions *workingConditions)
            : workingConditions(workingConditions), grid(setup.numerics.energyGrid), mixture(&grid),
              g_c(grid.cellNumber), eedf(grid.cellNumber), elasticMatrix(grid.cellNumber, grid.cellNumber),
              continuousMatrix(grid.cellNumber, grid.cellNumber) {

        mixture.initialize(setup, workingConditions);

        this->eedfType = setup.eedfType;
        this->shapeParameter = setup.shapeParameter;
        this->ionizationOperatorType = setup.ionizationOperatorType;
        this->growthModelType = setup.growthModelType;
        this->includeEECollisions = setup.includeEECollisions;

        this->plot("Total Elastic Cross Section N2", "Energy (eV)", "Cross Section (m^2)",
                   mixture.grid->getNodes(), mixture.totalCrossSection);
    }

    void ElectronKinetics::plot(const std::string &title, const std::string &xlabel, const std::string &ylabel,
                                const Vector &x, const Vector &y) {

        std::cout << "set term qt 1" << std::endl;
        std::cout << "unset key" << std::endl;
        std::cout << "set xlabel \"" << xlabel << "\"" << std::endl;
        std::cout << "set ylabel \"" << ylabel << "\"" << std::endl;
        std::cout << "set title \"" << title << "\"" << std::endl;
        std::cout << "set xrange [" << x[0] << ":" << x[x.size() - 1] << "]" << std::endl;
        std::cout << "plot '-' w l" << std::endl;
        for (uint32_t i = 0; i < x.size(); ++i) {
            std::cout << x[i] << "\t" << y[i] << '\n';
        }
        std::cout << "e" << std::endl;
    }
}