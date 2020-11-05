//
// Created by daan on 04-07-2019.
//

#include "LoKI-B/Output.h"
#include "LoKI-B/EedfCollision.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/StandardPaths.h"

#include <filesystem>
#include <fstream>
#include <iomanip>

namespace loki
{
namespace fs = std::filesystem;

Output::~Output()
{
}

Output::Output(const Setup &s, const WorkingConditions *workingConditions, const JobManager *jobManager)
    : workingConditions(workingConditions), jobManager(jobManager)
{

    for (const auto &entry : s.output.dataFiles)
    {
        if (entry == "eedf")
        {
            saveEedf = true;
        }
        else if (entry == "swarmParameters")
        {
            saveSwarm = true;
        }
        else if (entry == "rateCoefficients")
        {
            saveRates = true;
        }
        else if (entry == "powerBalance")
        {
            savePower = true;
        }
        else if (entry == "lookUpTable")
        {
            saveTable = true;
        }
    }
}

Output::Output(const json_type &cnf, const WorkingConditions *workingConditions, const JobManager *jobManager)
    : workingConditions(workingConditions),
      jobManager(jobManager)
{

    for (const auto &entry : cnf.at("output").at("dataFiles"))
    {
        if (entry == "eedf")
        {
            saveEedf = true;
        }
        else if (entry == "swarmParameters")
        {
            saveSwarm = true;
        }
        else if (entry == "rateCoefficients")
        {
            saveRates = true;
        }
        else if (entry == "powerBalance")
        {
            savePower = true;
        }
        else if (entry == "lookUpTable")
        {
            saveTable = true;
        }
    }
}

void Output::saveCycle(const Grid &energyGrid, const Vector &eedf, const WorkingConditions &wc, const Power &power,
                       const std::vector<EedfGas *> &gases, const SwarmParameters &swarmParameters,
                       const std::vector<RateCoefficient> &rateCoefficients,
                       const std::vector<RateCoefficient> &extraRateCoefficients, const Vector &firstAnisotropy)
{
    setDestination(jobManager->getCurrentJobFolder());
    if (saveEedf)
        writeEedf(eedf, firstAnisotropy, energyGrid.getCells());
    if (saveSwarm)
        writeSwarm(swarmParameters);
    if (savePower)
        writePower(power, gases);
    if (saveRates)
        writeRateCoefficients(rateCoefficients, extraRateCoefficients);
    if (saveTable)
        writeLookuptable(power, swarmParameters);
}

FileOutput::FileOutput(const Setup &setup, const WorkingConditions *workingConditions, const JobManager *jobManager)
 : Output(setup,workingConditions,jobManager), m_folder(OUTPUT "/" + setup.output.folder)
{
    createPath();
    std::ofstream ofs{m_folder + "/setup.in"};
    ofs << setup.fileContent << std::endl;
}

FileOutput::FileOutput(const json_type &cnf, const WorkingConditions *workingConditions, const JobManager *jobManager)
 : Output(cnf,workingConditions,jobManager), m_folder(OUTPUT "/" + cnf.at("output").at("folder").get<std::string>())
{
    createPath();
    std::ofstream ofs{m_folder + "/setup.json"};
    ofs << cnf.dump(1, '\t') << std::endl;
}

void FileOutput::setDestination(const std::string& subFolder)
{
    m_subFolder = subFolder;
    fs::path subPath(m_folder + '/' + m_subFolder);
    fs::create_directory(subPath);
}

void FileOutput::writeEedf(const Vector &eedf, const Vector &firstAnisotropy, const Vector &energies) const
{
    auto *file = std::fopen((m_folder + "/" + m_subFolder + "/eedf.txt").c_str(), "w");

    fprintf(file, "Energy (eV)          EEDF (eV^-(3/2))     First Anisotropy\n");

    for (uint32_t i = 0; i < energies.size(); ++i)
    {
        fprintf(file, "%.14e %.14e %.14e\n", energies[i], eedf[i], firstAnisotropy[i]);
    }

    fclose(file);
}

void FileOutput::writeSwarm(const SwarmParameters &swarmParameters) const
{
    auto *file = std::fopen((m_folder + "/" + m_subFolder + "/swarm_parameters.txt").c_str(), "w");

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

void FileOutput::writePower(const Power &power, const std::vector<EedfGas *> &gases) const
{
    auto *file = std::fopen((m_folder + "/" + m_subFolder + "/power_balance.txt").c_str(), "w");

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
    for (uint32_t i = 0; i < 73; ++i)
        fprintf(file, "-");
    fprintf(file, "\n");
    fprintf(file, "                       Power Balance = %#+.14e (eVm3/s)\n", power.balance);
    fprintf(file, "              Relative Power Balance = % #.14e%%\n\n", power.relativeBalance * 100);
    fprintf(file, "           Elastic collisions (gain) = %#+.14e (eVm3/s)\n", power.elasticGain);
    fprintf(file, "           Elastic collisions (loss) = %#+.14e (eVm3/s) +\n", power.elasticLoss);
    for (uint32_t i = 0; i < 73; ++i)
        fprintf(file, "-");
    fprintf(file, "\n");
    fprintf(file, "            Elastic collisions (net) = %#+.14e (eVm3/s)\n\n", power.elasticNet);
    fprintf(file, "                          CAR (gain) = %#+.14e (eVm3/s)\n", power.carGain);
    fprintf(file, "                          CAR (gain) = %#+.14e (eVm3/s) +\n", power.carLoss);
    for (uint32_t i = 0; i < 73; ++i)
        fprintf(file, "-");
    fprintf(file, "\n");
    fprintf(file, "                           CAR (net) = %#+.14e (eVm3/s)\n\n", power.carNet);
    fprintf(file, "     Excitation inelastic collisions = %#+.14e (eVm3/s)\n", power.excitationIne);
    fprintf(file, "  Excitation superelastic collisions = %#+.14e (eVm3/s) +\n", power.excitationSup);
    for (uint32_t i = 0; i < 73; ++i)
        fprintf(file, "-");
    fprintf(file, "\n");
    fprintf(file, "         Excitation collisions (net) = %#+.14e (eVm3/s)\n\n", power.excitationNet);
    fprintf(file, "    Vibrational inelastic collisions = %#+.14e (eVm3/s)\n", power.vibrationalIne);
    fprintf(file, " Vibrational superelastic collisions = %#+.14e (eVm3/s) +\n", power.vibrationalSup);
    for (uint32_t i = 0; i < 73; ++i)
        fprintf(file, "-");
    fprintf(file, "\n");
    fprintf(file, "        Vibrational collisions (net) = %#+.14e (eVm3/s)\n\n", power.vibrationalNet);
    fprintf(file, "     Rotational inelastic collisions = %#+.14e (eVm3/s)\n", power.rotationalIne);
    fprintf(file, "  Rotational superelastic collisions = %#+.14e (eVm3/s) +\n", power.rotationalSup);
    for (uint32_t i = 0; i < 73; ++i)
        fprintf(file, "-");
    fprintf(file, "\n");
    fprintf(file, "         Rotational collisions (net) = %#+.14e (eVm3/s)\n", power.rotationalNet);

    for (const auto &gas : gases)
    {
        const GasPower &gasPower = gas->getPower();

        fprintf(file, "\n");
        for (uint32_t i = 0; i < 37; ++i)
            fprintf(file, "*");
        fprintf(file, " %s ", gas->name.c_str());
        for (uint32_t i = 0; i < 39 - gas->name.length(); ++i)
            fprintf(file, "*");
        fprintf(file, "\n\n");

        fprintf(file, "     Excitation inelastic collisions = %#+.14e (eVm3/s)\n", gasPower.excitationIne);
        fprintf(file, "  Excitation superelastic collisions = %#+.14e (eVm3/s) +\n", gasPower.excitationSup);
        for (uint32_t i = 0; i < 73; ++i)
            fprintf(file, "-");
        fprintf(file, "\n");
        fprintf(file, "         Excitation collisions (net) = %#+.14e (eVm3/s)\n\n", gasPower.excitationNet);
        fprintf(file, "    Vibrational inelastic collisions = %#+.14e (eVm3/s)\n", gasPower.vibrationalIne);
        fprintf(file, " Vibrational superelastic collisions = %#+.14e (eVm3/s) +\n", gasPower.vibrationalSup);
        for (uint32_t i = 0; i < 73; ++i)
            fprintf(file, "-");
        fprintf(file, "\n");
        fprintf(file, "        Vibrational collisions (net) = %#+.14e (eVm3/s)\n\n", gasPower.vibrationalNet);
        fprintf(file, "     Rotational inelastic collisions = %#+.14e (eVm3/s)\n", gasPower.rotationalIne);
        fprintf(file, "  Rotational superelastic collisions = %#+.14e (eVm3/s) +\n", gasPower.rotationalSup);
        for (uint32_t i = 0; i < 73; ++i)
            fprintf(file, "-");
        fprintf(file, "\n");
        fprintf(file, "         Rotational collisions (net) = %#+.14e (eVm3/s)\n\n", gasPower.rotationalNet);
        fprintf(file, "               Ionization collisions = %#+.14e (eVm3/s)\n", gasPower.ionizationIne);
        fprintf(file, "               Attachment collisions = %#+.14e (eVm3/s)\n", gasPower.attachmentIne);
    }

    fclose(file);
}

void FileOutput::writeRateCoefficients(const std::vector<RateCoefficient> &rateCoefficients,
                                   const std::vector<RateCoefficient> &extraRateCoefficients) const
{
    auto *file = std::fopen((m_folder + "/" + m_subFolder + "/rate_coefficients.txt").c_str(), "w");

    fprintf(file, "Ine.R.Coeff.(m3/s)   Sup.R.Coeff.(m3/s)   Description\n");

    for (const auto &rateCoeff : rateCoefficients)
    {
        std::stringstream ss;

        ss << *rateCoeff.collision;

        fprintf(file, "%20.14e %20.14e %s\n", rateCoeff.inelastic, rateCoeff.superelastic, ss.rdbuf()->str().c_str());
    }

    fprintf(file, "\n");
    for (uint32_t i = 0; i < 27; ++i)
        fprintf(file, "-");
    fprintf(file, "\n* Extra Rate Coefficients *\n");
    for (uint32_t i = 0; i < 27; ++i)
        fprintf(file, "-");
    fprintf(file, "\n\n");

    for (const auto &rateCoeff : extraRateCoefficients)
    {
        std::stringstream ss;

        ss << *rateCoeff.collision;

        fprintf(file, "%20.14e %20.14e %s\n", rateCoeff.inelastic, rateCoeff.superelastic, ss.rdbuf()->str().c_str());
    }

    fclose(file);
}

void FileOutput::writeLookuptable(const Power &power, const SwarmParameters &swarmParameters) const
{
    std::FILE *file;

    if (initTable)
    {
        file = std::fopen((m_folder + "/lookup_table.txt").c_str(), "w");

        fprintf(file, "RedField(Td)         RedDif(1/(ms))       RedMob(1/(msV))      RedTow(m2)           "
                      "RedAtt(m2)           MeanE(eV)            CharE(eV)            EleTemp(eV)          "
                      "DriftVelocity(m/s)   RelativePowerBalance\n");

        fclose(file);
        initTable = false;
    }

    file = std::fopen((m_folder + "/lookup_table.txt").c_str(), "a");

    fprintf(file, "%20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %19.14e%%\n",
            workingConditions->reducedField, swarmParameters.redDiffCoeff, swarmParameters.redMobCoeff,
            swarmParameters.redTownsendCoeff, swarmParameters.redAttCoeff, swarmParameters.meanEnergy,
            swarmParameters.characEnergy, swarmParameters.Te, swarmParameters.driftVelocity,
            power.relativeBalance * 100);

    fclose(file);
}

void FileOutput::createPath()
{
    fs::path path(m_folder);

    if (fs::exists(path))
    {
        Log<Message>::Warning("The output folder \"" + m_folder + "\" already exists, results might be overriden.");

        std::string newFolder;

        simPathExists.emit(newFolder);

        if (!newFolder.empty())
        {
            m_folder = OUTPUT "/" + newFolder;
            path = fs::path(m_folder);

            fs::create_directories(path);
        }
    }
    else
    {
        fs::create_directories(path);
    }
}

JsonOutput::JsonOutput(json_type& root, const Setup &setup, const WorkingConditions *workingConditions, const JobManager *jobManager)
 : Output(setup,workingConditions,jobManager), m_root(root), m_active(nullptr)
{
    m_root["setup"] = setup.fileContent;
}

JsonOutput::JsonOutput(json_type& root, const json_type &cnf, const WorkingConditions *workingConditions, const JobManager *jobManager)
 : Output(cnf,workingConditions,jobManager), m_root(root), m_active(nullptr)
{
    m_root["setup"] = cnf;
}

void JsonOutput::setDestination(const std::string& subFolder)
{
    m_active = &m_root[subFolder];
}

void JsonOutput::writeEedf(const Vector &eedf, const Vector &firstAnisotropy, const Vector &energies) const
{
    json_type& out = (*m_active)["eedf"];
    out["labels"] = { "Energy", "EEDF", "First Anisotropy"};
    out["units"] = { "eV", "eV^-(3/2)", "1"};
    json_type& data = out["data"];
    for (uint32_t i = 0; i < energies.size(); ++i)
    {
        data.push_back(json_type{energies[i], eedf[i], firstAnisotropy[i]});
    }

}

void JsonOutput::writeSwarm(const SwarmParameters &swarmParameters) const
{
    json_type& out = (*m_active)["swarm_parameters"];
}

void JsonOutput::writePower(const Power &power, const std::vector<EedfGas *> &gases) const
{
    json_type& out = (*m_active)["power_balance"];
}

void JsonOutput::writeRateCoefficients(const std::vector<RateCoefficient> &rateCoefficients,
                   const std::vector<RateCoefficient> &extraRateCoefficients) const
{
    json_type& out = (*m_active)["rate_coefficients"];
}

void JsonOutput::writeLookuptable(const Power &power, const SwarmParameters &swarmParameters) const
{
    json_type& out = (*m_active)["lookup_table"];
}


} // namespace loki
