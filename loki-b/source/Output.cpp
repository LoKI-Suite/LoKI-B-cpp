//
// Created by daan on 04-07-2019.
//

#include "Output.h"
#include "Log.h"
#include "StandardPaths.h"

#include <filesystem>
#include <iomanip>

namespace loki {
    namespace fs = std::filesystem;

    Output::Output(const OutputSetup &s, const Grid *grid, const WorkingConditions *workingConditions,
                   const JobManager *jobManager)
            : folder(OUTPUT "/" + s.folder), grid(grid), workingConditions(workingConditions), jobManager(jobManager) {

        for (const auto &entry : s.dataFiles) {
            if (entry == "eedf") {
                saveEedf = true;
            } else if (entry == "swarmParameters") {
                saveSwarm = true;
            } else if (entry == "rateCoefficients") {
                saveRates = true;
            } else if (entry == "powerBalance") {
                savePower = true;
            } else if (entry == "lookUpTable") {
                saveTable = true;
            }
        }

        fs::path path(folder);

        if (fs::exists(path)) {
            simPathExists.emit();
        } else {
            fs::create_directories(path);
        }
    }

    void Output::saveCycle(const Vector &eedf, const Power &power, const std::vector<EedfGas *> &gasses,
                           const SwarmParameters &swarmParameters,
                           const std::vector<RateCoefficient> &rateCoefficients,
                           const std::vector<RateCoefficient> &extraRateCoefficients, const Vector &firstAnisotropy) {

        subFolder = jobManager->getCurrentJobFolder();

        fs::path subPath(folder + '/' + subFolder);
        fs::create_directory(subPath);

        if (saveEedf) writeEedf(eedf, firstAnisotropy, grid->getCells());
        if (saveSwarm) writeSwarm(swarmParameters);
        if (savePower) writePower(power, gasses);
        if (saveRates) writeRateCoefficients(rateCoefficients, extraRateCoefficients);
        if (saveTable) writeLookuptable(power, swarmParameters);

//        this->plot("Eedf", "Energy (eV)", "Eedf (Au)", grid->getCells(), eedf);
//        this->plot("First Anisotropy", "Energy (eV)", "First Anisotropy (Au)", grid.getCells(), firstAnisotropy);
    }

    void Output::writeEedf(const Vector &eedf, const Vector &firstAnisotropy, const Vector &energies) {
        auto *file = std::fopen((folder + "/" + subFolder + "/eedf.txt").c_str(), "w");

        fprintf(file, "Energy (eV)          EEDF (eV^-(3/2))     First Anisotropy\n");

        for (uint32_t i = 0; i < grid->cellNumber; ++i) {
            fprintf(file, "%.14e %.14e %.14e\n", energies[i], eedf[i], firstAnisotropy[i]);
        }

        fclose(file);
    }

    void Output::writeSwarm(const SwarmParameters &swarmParameters) {
        auto *file = std::fopen((folder + "/" + subFolder + "/swarm_parameters.txt").c_str(), "w");

        fprintf(file, "         Reduced electric field = %#.14e (Td)\n", workingConditions->reducedField);
        fprintf(file, "  Reduced diffusion coefficient = %#.14e (1/(ms))\n", swarmParameters.redDiffCoeff);
        fprintf(file, "   Reduced mobility coefficient = %#.14e (1/(msV))\n", swarmParameters.redMobCoeff);
        fprintf(file, "   Reduced Townsend coefficient = %#.14e (m2)\n", swarmParameters.redTownsendCoeff);
        fprintf(file, " Reduced attachment coefficient = %#.14e (m2)\n", swarmParameters.redAttCoeff);
        fprintf(file, "                    Mean energy = %#.14e (eV)\n", swarmParameters.meanEnergy);
        fprintf(file, "          Characteristic energy = %#.14e (eV)\n", swarmParameters.characEnergy);
        fprintf(file, "           Electron temperature = %#.14e (eV)\n", swarmParameters.Te);
        fprintf(file, "                 Drift velocity = %#.14e (m/s)\n", swarmParameters.driftVelocity);

        fclose(file);
    }

    void Output::writePower(const Power &power, const std::vector<EedfGas *> &gasses) {
        auto *file = std::fopen((folder + "/" + subFolder + "/power_balance.txt").c_str(), "w");

        fprintf(file, "                               Field = %#+.14e (eVm3/s)\n", power.field);
        fprintf(file, "           Elastic collisions (gain) = %#+.14e (eVm3/s)\n", power.elasticGain);
        fprintf(file, "           Elastic collisions (loss) = %#+.14e (eVm3/s)\n", power.elasticLoss);
        fprintf(file, "                          CAR (gain) = %#+.14e (eVm3/s)\n", power.carGain);
        fprintf(file, "                          CAR (loss) = %#+.14e (eVm3/s)\n", power.carLoss);
        fprintf(file, "     Excitation inelastic collisions = %#+.14e (eVm3/s)\n", power.excitationIne);
        fprintf(file, "  Excitation superelastic collisions = %#+.14e (eVm3/s)\n", power.excitationSup);
        fprintf(file, "    Vibrational inelastic collisions = %#+.14e (eVm3/s)\n", power.vibrationalIne);
        fprintf(file, " Vibrational superelastic collisions = %#+.14e (eVm3/s)\n", power.vibrationalSup);
        fprintf(file, "     Rotational inelastic collisions = %#+.14e (eVm3/s)\n", power.rotationalIne);
        fprintf(file, "  Rotational superelastic collisions = %#+.14e (eVm3/s)\n", power.rotationalSup);
        fprintf(file, "               Ionization collisions = %#+.14e (eVm3/s)\n", power.ionizationIne);
        fprintf(file, "               Attachment collisions = %#+.14e (eVm3/s)\n", power.attachmentIne);
        fprintf(file, "             Electron density growth = %#+.14e (eVm3/s) +\n", power.eDensGrowth);
        for (uint32_t i = 0; i < 73; ++i) fprintf(file, "-");
        fprintf(file, "\n");
        fprintf(file, "                       Power Balance = %#+.14e (eVm3/s)\n", power.balance);
        fprintf(file, "              Relative Power Balance = % #.14e%%\n\n", power.relativeBalance * 100);
        fprintf(file, "           Elastic collisions (gain) = %#+.14e (eVm3/s)\n", power.elasticGain);
        fprintf(file, "           Elastic collisions (loss) = %#+.14e (eVm3/s) +\n", power.elasticLoss);
        for (uint32_t i = 0; i < 73; ++i) fprintf(file, "-");
        fprintf(file, "\n");
        fprintf(file, "            Elastic collisions (net) = %#+.14e (eVm3/s)\n\n", power.elasticNet);
        fprintf(file, "                          CAR (gain) = %#+.14e (eVm3/s)\n", power.carGain);
        fprintf(file, "                          CAR (gain) = %#+.14e (eVm3/s) +\n", power.carLoss);
        for (uint32_t i = 0; i < 73; ++i) fprintf(file, "-");
        fprintf(file, "\n");
        fprintf(file, "                           CAR (net) = %#+.14e (eVm3/s)\n\n", power.carNet);
        fprintf(file, "     Excitation inelastic collisions = %#+.14e (eVm3/s)\n", power.excitationIne);
        fprintf(file, "  Excitation superelastic collisions = %#+.14e (eVm3/s) +\n", power.excitationSup);
        for (uint32_t i = 0; i < 73; ++i) fprintf(file, "-");
        fprintf(file, "\n");
        fprintf(file, "         Excitation collisions (net) = %#+.14e (eVm3/s)\n\n", power.excitationNet);
        fprintf(file, "    Vibrational inelastic collisions = %#+.14e (eVm3/s)\n", power.vibrationalIne);
        fprintf(file, " Vibrational superelastic collisions = %#+.14e (eVm3/s) +\n", power.vibrationalSup);
        for (uint32_t i = 0; i < 73; ++i) fprintf(file, "-");
        fprintf(file, "\n");
        fprintf(file, "        Vibrational collisions (net) = %#+.14e (eVm3/s)\n\n", power.vibrationalNet);
        fprintf(file, "     Rotational inelastic collisions = %#+.14e (eVm3/s)\n", power.rotationalIne);
        fprintf(file, "  Rotational superelastic collisions = %#+.14e (eVm3/s) +\n", power.rotationalSup);
        for (uint32_t i = 0; i < 73; ++i) fprintf(file, "-");
        fprintf(file, "\n");
        fprintf(file, "         Rotational collisions (net) = %#+.14e (eVm3/s)\n", power.rotationalNet);

        for (const auto *gas : gasses) {
            const GasPower &gasPower = gas->getPower();

            fprintf(file, "\n");
            for (uint32_t i = 0; i < 37; ++i) fprintf(file, "*");
            fprintf(file, " %s ", gas->name.c_str());
            for (uint32_t i = 0; i < 39 - gas->name.length(); ++i) fprintf(file, "*");
            fprintf(file, "\n\n");

            fprintf(file, "     Excitation inelastic collisions = %#+.14e (eVm3/s)\n", gasPower.excitationIne);
            fprintf(file, "  Excitation superelastic collisions = %#+.14e (eVm3/s) +\n", gasPower.excitationSup);
            for (uint32_t i = 0; i < 73; ++i) fprintf(file, "-");
            fprintf(file, "\n");
            fprintf(file, "         Excitation collisions (net) = %#+.14e (eVm3/s)\n\n", gasPower.excitationNet);
            fprintf(file, "    Vibrational inelastic collisions = %#+.14e (eVm3/s)\n", gasPower.vibrationalIne);
            fprintf(file, " Vibrational superelastic collisions = %#+.14e (eVm3/s) +\n", gasPower.vibrationalSup);
            for (uint32_t i = 0; i < 73; ++i) fprintf(file, "-");
            fprintf(file, "\n");
            fprintf(file, "        Vibrational collisions (net) = %#+.14e (eVm3/s)\n\n", gasPower.vibrationalNet);
            fprintf(file, "     Rotational inelastic collisions = %#+.14e (eVm3/s)\n", gasPower.rotationalIne);
            fprintf(file, "  Rotational superelastic collisions = %#+.14e (eVm3/s) +\n", gasPower.rotationalSup);
            for (uint32_t i = 0; i < 73; ++i) fprintf(file, "-");
            fprintf(file, "\n");
            fprintf(file, "         Rotational collisions (net) = %#+.14e (eVm3/s)\n\n", gasPower.rotationalNet);
            fprintf(file, "               Ionization collisions = %#+.14e (eVm3/s)\n", gasPower.ionizationIne);
            fprintf(file, "               Attachment collisions = %#+.14e (eVm3/s)\n", gasPower.attachmentIne);
        }

        fclose(file);
    }

    void Output::writeRateCoefficients(const std::vector<RateCoefficient> &rateCoefficients,
                                       const std::vector<RateCoefficient> &extraRateCoefficients) {
        auto *file = std::fopen((folder + "/" + subFolder + "/rate_coefficients.txt").c_str(), "w");

        fprintf(file, "Ine.R.Coeff.(m3/s)   Sup.R.Coeff.(m3/s)   Description\n");

        for (const auto &rateCoeff : rateCoefficients) {
            std::stringstream ss;

            ss << *rateCoeff.collision;

            fprintf(file, "%20.14e %20.14e %s\n", rateCoeff.inelastic, rateCoeff.superelastic,
                    ss.rdbuf()->str().c_str());
        }

        fprintf(file, "\n");
        for (uint32_t i = 0; i < 27; ++i) fprintf(file, "-");
        fprintf(file, "\n* Extra Rate Coefficients *\n");
        for (uint32_t i = 0; i < 27; ++i) fprintf(file, "-");
        fprintf(file, "\n\n");

        for (const auto &rateCoeff : extraRateCoefficients) {
            std::stringstream ss;

            ss << *rateCoeff.collision;

            fprintf(file, "%20.14e %20.14e %s\n", rateCoeff.inelastic, rateCoeff.superelastic,
                    ss.rdbuf()->str().c_str());
        }

        fclose(file);
    }

    void Output::writeLookuptable(const Power &power, const SwarmParameters &swarmParameters) {
        std::FILE *file;

        if (initTable) {
            file = std::fopen((folder + "/lookup_table.txt").c_str(), "w");

            fprintf(file, "RedField(Td)         RedDif(1/(ms))       RedMob(1/(msV))      RedTow(m2)           "
                          "RedAtt(m2)           MeanE(eV)            CharE(eV)            EleTemp(eV)          "
                          "DriftVelocity(m/s)   RelativePowerBalance\n");

            fclose(file);
            initTable = false;
        }

        file = std::fopen((folder + "/lookup_table.txt").c_str(), "a");

        fprintf(file, "%20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %19.14e%%\n",
                workingConditions->reducedField, swarmParameters.redDiffCoeff, swarmParameters.redMobCoeff,
                swarmParameters.redTownsendCoeff, swarmParameters.redAttCoeff, swarmParameters.meanEnergy,
                swarmParameters.characEnergy, swarmParameters.Te, swarmParameters.driftVelocity,
                power.relativeBalance * 100);

        fclose(file);
    }

    void Output::plot(const std::string &title, const std::string &xlabel, const std::string &ylabel,
                      const Vector &x, const Vector &y) {
        std::cout << "unset key" << std::endl;
        std::cout << "set xlabel \"" << xlabel << "\"" << std::endl;
        std::cout << "set ylabel \"" << ylabel << "\"" << std::endl;
        std::cout << "set title \"" << title << "\"" << std::endl;
        std::cout << "set xrange [" << x[0] << ":" << x[x.size() - 1] << "]" << std::endl;
        std::cout << "set logscale y" << std::endl;
        std::cout << "plot '-' w l" << std::endl;
        for (uint32_t i = 0; i < x.size(); ++i) {
            std::cout << x[i] << "\t" << y[i] << '\n';
        }
        std::cout << "e" << std::endl;
    }


}