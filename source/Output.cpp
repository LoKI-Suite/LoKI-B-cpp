/** \file
 *
 *  Implementations of classes for output generation.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2024 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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
 *  \author Jan van Dijk, Daan Boer and Jop Hendrikx
 *  \date   4 July 2019
 */

#include "LoKI-B/Output.h"
#include "LoKI-B/EedfCollisions.h"
#include "LoKI-B/EedfMixture.h"
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

Output::Output(const json_type &cnf, const WorkingConditions *workingConditions)
    : m_workingConditions(workingConditions),
    m_isBoltzmann(cnf.at("electronKinetics").at("eedfType")=="boltzmann"),
    m_isSimulationHF(m_workingConditions->reducedExcFreqSI()>0.0)
{
    saveEedf = false;
    savePower = false;
    saveSwarm = false;
    saveRates = false;
    saveTable = false;

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
                       const EedfCollisionDataMixture &collData, const SwarmParameters &swarmParameters,
                       const Vector *firstAnisotropy)
{
    setDestination(wc.getCurrentJobFolder());
    if (saveEedf)
        writeEedf(eedf, firstAnisotropy, energyGrid.getCells());
    if (saveSwarm)
        writeSwarm(swarmParameters);
    if (savePower)
        writePower(power, collData);
    if (saveRates)
        writeRateCoefficients(collData.m_rateCoefficients, collData.m_rateCoefficientsExtra);
    if (saveTable)
        writeLookuptable(power, swarmParameters);
}

FileOutput::FileOutput(const json_type &cnf, const WorkingConditions *workingConditions,
        const PathExistsHandler& handler)
 : Output(cnf,workingConditions), m_folder(LOKIB_OUTPUT_DIR "/" + cnf.at("output").at("folder").get<std::string>())
{
    m_initTable = true;
    createPath(handler);
    std::ofstream ofs{m_folder + "/setup.json"};
    ofs << cnf.dump(1, '\t') << std::endl;
}

void FileOutput::setDestination(const std::string& subFolder)
{
    m_subFolder = subFolder;
    fs::path subPath(m_folder + '/' + m_subFolder);
    fs::create_directory(subPath);
}

void FileOutput::writeTerm(std::ostream& os, const std::string& name, const std::string& unit, double value, bool plus) const
{
    os << std::setw(36) << name << " = "
        << std::showpos << std::scientific << std::setprecision(14) << value
        << " (" << unit << ")";
    if (plus) os << " +";
    os << std::endl;
}

void FileOutput::writeEedf(const Vector &eedf, const Vector *firstAnisotropy, const Vector &energies) const
{
    std::ofstream os(m_folder + "/" + m_subFolder + "/eedf.txt");
    os << "Energy (eV)          EEDF (eV^-(3/2))";
    if (firstAnisotropy)
        os << "     First Anisotropy" << std::endl;
    for (Vector::Index i = 0; i < energies.size(); ++i)
    {
        os << std::scientific << std::setprecision(14) << energies[i];
        os << ' ';
        os << std::scientific << std::setprecision(14) << eedf[i];
        if (firstAnisotropy)
        {
            os << ' ';
            os << std::scientific << std::setprecision(14) << (*firstAnisotropy)[i];
        }
        os << std::endl;
    }
}

void FileOutput::writeSwarm(const SwarmParameters &swarmParameters) const
{
    std::ofstream os(m_folder + "/" + m_subFolder + "/swarm_parameters.txt");
    writeTerm(os,"Reduced electric field", "Td", m_workingConditions->reducedField());
    writeTerm(os,"Reduced diffusion coefficient","1/(ms)",swarmParameters.redDiffCoeff);
    writeTerm(os,"Reduced mobility coefficient","1/(msV)",swarmParameters.redMobCoeff);
    if (isSimulationHF())
    {
        writeTerm(os,"Re(Reduced mobility HF)","1/(msV)", swarmParameters.redMobilityHF.real());
        writeTerm(os,"Im(Reduced mobility HF)","1/(msV)", swarmParameters.redMobilityHF.imag());
    }
    else
    {
        writeTerm(os,"Drift velocity","m/s",swarmParameters.driftVelocity);
        writeTerm(os,"Reduced Townsend coefficient","m2",swarmParameters.redTownsendCoeff);
        writeTerm(os,"Reduced attachment coefficient","m2",swarmParameters.redAttCoeff);
    }
    writeTerm(os,"Reduced energy diffusion coefficient","eV/(ms)",swarmParameters.redDiffCoeffEnergy);
    writeTerm(os,"Reduced energy mobility","eV/(msV)",swarmParameters.redMobilityEnergy);
    writeTerm(os,"Mean energy","eV",swarmParameters.meanEnergy);
    writeTerm(os,"Characteristic energy","eV",swarmParameters.characEnergy);
    writeTerm(os,"Electron temperature","eV",swarmParameters.Te);
}

void FileOutput::writePower(const Power &power, const EedfCollisionDataMixture& collData) const
{
    std::ofstream os(m_folder + "/" + m_subFolder + "/power_balance.txt");
    writeTerm(os, "Field","eVm3/s", power.field);
    writeTerm(os, "Elastic collisions (gain)","eVm3/s", power.elasticGain);
    writeTerm(os, "Elastic collisions (loss)","eVm3/s", power.elasticLoss);
    writeTerm(os, "CAR (gain)","eVm3/s", power.carGain);
    writeTerm(os, "CAR (loss)","eVm3/s", power.carLoss);
    writeTerm(os, "Excitation inelastic collisions","eVm3/s", power.excitation.forward);
    writeTerm(os, "Excitation superelastic collisions","eVm3/s", power.excitation.backward);
    writeTerm(os, "Vibrational inelastic collisions","eVm3/s", power.vibrational.forward);
    writeTerm(os, "Vibrational superelastic collisions","eVm3/s", power.vibrational.backward);
    writeTerm(os, "Rotational inelastic collisions","eVm3/s", power.rotational.forward);
    writeTerm(os, "Rotational superelastic collisions","eVm3/s", power.rotational.backward);
    writeTerm(os, "Ionization collisions","eVm3/s", power.ionization.forward); // no recombination
    writeTerm(os, "Attachment collisions","eVm3/s", power.attachment.forward); // no detachment
    writeTerm(os, "Electron density growth","eVm3/s", power.eDensGrowth,true);
    os << std::string(73,'-') << std::endl;

    writeTerm(os,"Power Balance","eVm3/s", power.balance);
    writeTerm(os,"Relative Power Balance", "%", power.relativeBalance * 100);
    writeTerm(os,"Elastic collisions (gain)","eVm3/s", power.elasticGain);
    writeTerm(os,"Elastic collisions (loss)","eVm3/s", power.elasticLoss,true);
    writeTerm(os,"Elastic electron-electron","eVm3/s", power.electronElectron,true);
    os << std::string(73,'-') << std::endl;
    writeTerm(os,"Elastic collisions (net)","eVm3/s", power.elasticNet);
    writeTerm(os,"CAR (gain)","eVm3/s", power.carGain);
    writeTerm(os,"CAR (loss)","eVm3/s", power.carLoss,true);
    os << std::string(73,'-') << std::endl;
    writeTerm(os,"CAR (net)","eVm3/s", power.carNet);
    writeTerm(os,"Excitation inelastic collisions","eVm3/s", power.excitation.forward);
    writeTerm(os,"Excitation superelastic collisions","eVm3/s", power.excitation.backward,true);
    os << std::string(73,'-') << std::endl;
    writeTerm(os,"Excitation collisions (net)","eVm3/s", power.excitation.net());
    writeTerm(os,"Vibrational inelastic collisions","eVm3/s", power.vibrational.forward);
    writeTerm(os,"Vibrational superelastic collisions","eVm3/s", power.vibrational.backward,true);
    os << std::string(73,'-') << std::endl;
    writeTerm(os,"Vibrational collisions (net)","eVm3/s", power.vibrational.net());
    writeTerm(os,"Rotational inelastic collisions","eVm3/s", power.rotational.forward);
    writeTerm(os,"Rotational superelastic collisions","eVm3/s", power.rotational.backward,true);
    os << std::string(73,'-') << std::endl;
    writeTerm(os,"Rotational collisions (net)","eVm3/s", power.rotational.net());

    for (const auto &cd : collData.data_per_gas())
    {
        const GasPower &gasPower = cd.getPower();

        os << std::endl;
        os << std::string(37,'*') << ' ' << cd.gas().name() << ' ' << std::string(39 - cd.gas().name().length(),'*') << std::endl;
        os << std::endl;

        writeTerm(os,"Excitation inelastic collisions","eVm3/s", gasPower.excitation.forward);
        writeTerm(os,"Excitation superelastic collisions","eVm3/s", gasPower.excitation.backward,true);
        os << std::string(73,'-') << std::endl;
        writeTerm(os,"Excitation collisions (net)","eVm3/s", gasPower.excitation.net());
        writeTerm(os,"Vibrational inelastic collisions","eVm3/s", gasPower.vibrational.forward);
        writeTerm(os,"Vibrational superelastic collisions","eVm3/s", gasPower.vibrational.backward,true);
        os << std::string(73,'-') << std::endl;
        writeTerm(os,"Vibrational collisions (net)","eVm3/s", gasPower.vibrational.net());
        writeTerm(os,"Rotational inelastic collisions","eVm3/s", gasPower.rotational.forward);
        writeTerm(os,"Rotational superelastic collisions","eVm3/s", gasPower.rotational.backward,true);
        os << std::string(73,'-') << std::endl;
        writeTerm(os,"Rotational collisions (net)","eVm3/s", gasPower.rotational.net());
        writeTerm(os,"Ionization collisions","eVm3/s", gasPower.ionization.forward); // no recombination
        writeTerm(os,"Attachment collisions","eVm3/s", gasPower.attachment.forward); // no detachment
    }
}

void FileOutput::writeRateCoefficients(const std::vector<RateCoefficient> &rateCoefficients,
                                   const std::vector<RateCoefficient> &extraRateCoefficients) const
{
    std::ofstream os(m_folder + "/" + m_subFolder + "/rate_coefficients.txt");
    os << "Ine.R.Coeff.(m3/s)   Sup.R.Coeff.(m3/s)   Description" << std::endl;
    for (const auto &rateCoeff : rateCoefficients)
    {
        os << std::setw(20) << std::scientific << std::setprecision(14) << rateCoeff.inelastic;
        os << ' ';
        os << std::setw(20) << std::scientific << std::setprecision(14) << rateCoeff.superelastic;
        os << ' ';
        os << *rateCoeff.collision;
        os << std::endl;
    }
    os << std::endl;
    os << std::string(27,'-');
    os << std::endl;
    os << "* Extra Rate Coefficients *" << std::endl;
    os << std::string(27,'-');
    os << std::endl;
    os << std::endl;
    for (const auto &rateCoeff : extraRateCoefficients)
    {
        os << std::setw(20) << std::scientific << std::setprecision(14) << rateCoeff.inelastic;
        os << ' ';
        os << std::setw(20) << std::scientific << std::setprecision(14) << rateCoeff.superelastic;
        os << ' ';
        os << *rateCoeff.collision;
        os << std::endl;
    }
}

void FileOutput::writeLookuptable(const Power &power, const SwarmParameters &swarmParameters) const
{
    /// \todo Write lookUpTablePower.txt and lookUpTableRateCoeff.txt
    std::ofstream os(m_folder + "/lookup_table.txt", m_initTable ? std::ios_base::trunc : std::ios_base::app);
    if (isBoltzmann())
    {
        if (m_initTable)
        {
            os
                << "RedField(Td)          "
                << "RedDif(1/(ms))        "
                << "RedMob(1/(msV))       ";
            if (isSimulationHF())
            {
                os
                    << "R[RedMobHF]((msV)^-1) "
                    << "I[RedMobHF]((msV)^-1)";
            }
            else
            {
                os
                    << "DriftVelocity(m/s)    "
                    << "RedTow(m2)            "
                    << "RedAtt(m2)            ";
            }
            os
                << "RedDiffE(eV/(ms))     "
                << "RedMobE(eV/(msV))     "
                << "MeanE(eV)             "
                << "CharE(eV)             "
                << "EleTemp(eV)           "
                << "RelativePowerBalance"
                << std::endl;
           m_initTable = false;
        }
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << m_workingConditions->reducedField();
        os << ' ';
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.redDiffCoeff;
        os << ' ';
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.redMobCoeff;
        if (isSimulationHF())
        {
            os << ' ';
            os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.redMobilityHF.real();
            os << ' ';
            os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.redMobilityHF.imag();
        }
        else
        {
            os << ' ';
            os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.driftVelocity;
            os << ' ';
            os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.redTownsendCoeff;
            os << ' ';
            os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.redAttCoeff;
        }
        os << ' ';
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.redDiffCoeffEnergy;
        os << ' ';
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.redMobilityEnergy;
        os << ' ';
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.meanEnergy;
        os << ' ';
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.characEnergy;
        os << ' ';
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.Te;
        os << ' ';
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << (power.relativeBalance * 100) << '%';
        os << std::endl;
    }
    else // prescribed EEDF
    {
        if (m_initTable)
        {
            os
                << "EleTemp(eV)           "
                << "RedField(Td)          "
                << "RedDif(1/(ms))        "
                << "RedMob(1/(msV))       ";
            if (isSimulationHF())
            {
                os
                    << "R[RedMobHF]((msV)^-1) "
                    << "I[RedMobHF]((msV)^-1)";
            }
            else
            {
                os
                    << "DriftVelocity(m/s)    "
                    << "RedTow(m2)            "
                    << "RedAtt(m2)            ";
            }
            os
                << "RedDiffE(eV/(ms))     "
                << "RedMobE(eV/(msV))     "
                << "MeanE(eV)             "
                << "CharE(eV)             "
                << "RelativePowerBalance"
                << std::endl;
           m_initTable = false;
        }
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.Te;
        os << ' ';
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << m_workingConditions->reducedField();
        os << ' ';
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.redDiffCoeff;
        os << ' ';
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.redMobCoeff;
        if (isSimulationHF())
        {
            os << ' ';
            os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.redMobilityHF.real();
            os << ' ';
            os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.redMobilityHF.imag();
        }
        else
        {
            os << ' ';
            os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.driftVelocity;
            os << ' ';
            os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.redTownsendCoeff;
            os << ' ';
            os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.redAttCoeff;
        }
        os << ' ';
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.redDiffCoeffEnergy;
        os << ' ';
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.redMobilityEnergy;
        os << ' ';
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.meanEnergy;
        os << ' ';
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << swarmParameters.characEnergy;
        os << ' ';
        os << std::showpos << std::setw(20) << std::scientific << std::setprecision(14) << (power.relativeBalance * 100) << '%';
        os << std::endl;
    }
}

void FileOutput::createPath(const PathExistsHandler& handler)
{
    fs::path path(m_folder);

    if (fs::exists(path))
    {
        Log<Message>::Warning("The output folder \"" + m_folder + "\" already exists, results might be overriden.");

        std::string newFolder;

        handler(newFolder);

        if (!newFolder.empty())
        {
            m_folder = LOKIB_OUTPUT_DIR "/" + newFolder;
            path = fs::path(m_folder);

            fs::create_directories(path);
        }
    }
    else
    {
        fs::create_directories(path);
    }
}

JsonOutput::JsonOutput(json_type& root, const json_type &cnf, const WorkingConditions *workingConditions)
 : Output(cnf,workingConditions), m_root(root), m_active(nullptr)
{
    m_root["setup"] = cnf;
}

void JsonOutput::setDestination(const std::string& subFolder)
{
    m_active = &m_root[subFolder];
}

json_type JsonOutput::makeQuantity(const std::string& name, double value, const std::string unit)
{
    return { { name, { { "value", value }, { "unit", unit } } } };
}

void JsonOutput::writeEedf(const Vector &eedf, const Vector *firstAnisotropy, const Vector &energies) const
{
    json_type& out = (*m_active)["eedf"];
    if (firstAnisotropy)
    {
        out["labels"] = { "Energy", "EEDF", "First Anisotropy"};
        out["units"] = { "eV", "eV^-(3/2)", "1"};
    }
    else
    {
        out["labels"] = { "Energy", "EEDF"};
        out["units"] = { "eV", "eV^-(3/2)"};
    }
    json_type& data = out["data"];
    for (Vector::Index i = 0; i < energies.size(); ++i)
    {
        if (firstAnisotropy)
        {
            data.push_back(json_type{energies[i], eedf[i], (*firstAnisotropy)[i]});
        }
        else
        {
            data.push_back(json_type{energies[i], eedf[i]});
        }
    }
}

void JsonOutput::writeSwarm(const SwarmParameters &swarmParameters) const
{
    /// \todo Handle isSimulationHF()
    json_type& out = (*m_active)["swarm_parameters"];
    out.push_back( makeQuantity("Reduced electric field", m_workingConditions->reducedField(), "Td") );
    out.push_back( makeQuantity("Reduced diffusion coefficient", swarmParameters.redDiffCoeff, "1/(m*s)") );
    out.push_back( makeQuantity("Reduced mobility coefficient", swarmParameters.redMobCoeff, "1/(m*s*V)") );
    out.push_back( makeQuantity("Reduced Townsend coefficient", swarmParameters.redTownsendCoeff, "m^2") );
    out.push_back( makeQuantity("Reduced attachment coefficient", swarmParameters.redAttCoeff, "m^2") );
    out.push_back( makeQuantity("Mean energy", swarmParameters.meanEnergy, "eV") );
    out.push_back( makeQuantity("Characteristic energy", swarmParameters.characEnergy, "eV") );
    out.push_back( makeQuantity("Electron temperature", swarmParameters.Te, "eV") );
    out.push_back( makeQuantity("Drift velocity", swarmParameters.driftVelocity, "m/s") );
}

void JsonOutput::writePower(const Power &power, const EedfCollisionDataMixture& collData) const
{
    /// \todo Handle isSimulationHF()
    json_type& out = (*m_active)["power_balance"];
    out.push_back( makeQuantity("Field", power.field, "eV*m^3/s") );
    //out.push_back( { "Field", "eV*m^3/s", power.field });
    out.push_back( makeQuantity("Elastic collisions (gain)", power.elasticGain, "eV*m^3/s") );
    out.push_back( makeQuantity("Elastic collisions (loss)", power.elasticLoss, "eV*m^3/s") );
    out.push_back( makeQuantity("CAR (gain)", power.carGain, "eV*m^3/s") );
    out.push_back( makeQuantity("CAR (loss)", power.carLoss, "eV*m^3/s") );
    out.push_back( makeQuantity("Excitation inelastic collisions", power.excitation.forward, "eV*m^3/s") );
    out.push_back( makeQuantity("Excitation superelastic collisions", power.excitation.backward, "eV*m^3/s") );
    out.push_back( makeQuantity("Vibrational inelastic collisions", power.vibrational.forward, "eV*m^3/s") );
    out.push_back( makeQuantity("Vibrational superelastic collisions", power.vibrational.backward, "eV*m^3/s") );
    out.push_back( makeQuantity("Rotational inelastic collisions", power.rotational.forward, "eV*m^3/s") );
    out.push_back( makeQuantity("Rotational superelastic collisions", power.rotational.backward, "eV*m^3/s") );
    out.push_back( makeQuantity("Ionization collisions", power.ionization.forward, "eV*m^3/s") ); // no recombination
    out.push_back( makeQuantity("Attachment collisions", power.attachment.forward, "eV*m^3/s") ); // no detachment
    out.push_back( makeQuantity("Electron density growth", power.eDensGrowth, "eV*m^3/s") );

    out.push_back( makeQuantity("Power Balance", power.balance, "eV*m^3/s") );
    out.push_back( makeQuantity("Relative Power Balance", power.relativeBalance * 100, "%") );
    out.push_back( makeQuantity("Elastic collisions (gain)", power.elasticGain, "eV*m^3/s") );
    out.push_back( makeQuantity("Elastic collisions (loss)", power.elasticLoss, "eV*m^3/s") );

    out.push_back( makeQuantity("Elastic collisions (net)", power.elasticNet, "eV*m^3/s") );
    out.push_back( makeQuantity("CAR (gain)", power.carGain, "eV*m^3/s") );
    out.push_back( makeQuantity("CAR (loss)", power.carLoss, "eV*m^3/s") );

    out.push_back( makeQuantity("CAR (net)", power.carNet, "eV*m^3/s") );
    out.push_back( makeQuantity("Excitation inelastic collisions", power.excitation.forward, "eV*m^3/s") );
    out.push_back( makeQuantity("Excitation superelastic collisions", power.excitation.backward, "eV*m^3/s") );

    out.push_back( makeQuantity("Excitation collisions (net)", power.excitation.net(), "eV*m^3/s") );
    out.push_back( makeQuantity("Vibrational inelastic collisions", power.vibrational.forward, "eV*m^3/s") );
    out.push_back( makeQuantity("Vibrational superelastic collisions", power.vibrational.backward, "eV*m^3/s") );

    out.push_back( makeQuantity("Vibrational collisions (net)", power.vibrational.net(), "eV*m^3/s") );
    out.push_back( makeQuantity("Rotational inelastic collisions", power.rotational.forward, "eV*m^3/s") );
    out.push_back( makeQuantity("Rotational superelastic collisions", power.rotational.backward, "eV*m^3/s") );

    out.push_back( makeQuantity("Rotational collisions (net)", power.rotational.net(), "eV*m^3/s") );

    for (const auto &cd : collData.data_per_gas())
    {
        const GasPower &gasPower = cd.getPower();
	json_type gas_out;
        gas_out.push_back( { { "name", cd.gas().name() } } );
        gas_out.push_back( makeQuantity("Excitation inelastic collisions", gasPower.excitation.forward, "eV*m^3/s") );
        gas_out.push_back( makeQuantity("Excitation superelastic collisions", gasPower.excitation.backward, "eV*m^3/s") );
        gas_out.push_back( makeQuantity("Excitation collisions (net)", gasPower.excitation.net(), "eV*m^3/s") );
        gas_out.push_back( makeQuantity("Vibrational inelastic collisions", gasPower.vibrational.forward, "eV*m^3/s") );
        gas_out.push_back( makeQuantity("Vibrational superelastic collisions", gasPower.vibrational.backward, "eV*m^3/s") );
        gas_out.push_back( makeQuantity("Vibrational collisions (net)", gasPower.vibrational.net(), "eV*m^3/s") );
        gas_out.push_back( makeQuantity("Rotational inelastic collisions", gasPower.rotational.forward, "eV*m^3/s") );
        gas_out.push_back( makeQuantity("Rotational superelastic collisions", gasPower.rotational.backward, "eV*m^3/s") );
        gas_out.push_back( makeQuantity("Rotational collisions (net)", gasPower.rotational.net(), "eV*m^3/s") );
        gas_out.push_back( makeQuantity("Ionization collisions", gasPower.ionization.forward, "eV*m^3/s") ); // no recombination
        gas_out.push_back( makeQuantity("Attachment collisions", gasPower.attachment.forward, "eV*m^3/s") ); // no detachment

        out.push_back( { "gas", gas_out } );
    }
}

void JsonOutput::writeRateCoefficients(const std::vector<RateCoefficient> &rateCoefficients,
                   const std::vector<RateCoefficient> &extraRateCoefficients) const
{
    /// \todo Handle isSimulationHF()
    if (rateCoefficients.size())
    {
        json_type& out = (*m_active)["rate_coefficients"];
        out["labels"] = { "Ine.R.Coeff.", "Sup.R.Coeff.", "Description"};
        out["units"] = { "m^3/s", "m^3/s", ""};
        json_type& data = out["data"];
        for (const auto &rateCoeff : rateCoefficients)
        {
            std::stringstream ss;
            ss << *rateCoeff.collision;
            data.push_back( json_type{rateCoeff.inelastic, rateCoeff.superelastic, ss.str() } );
        }
    }
    /** \todo See if we can merge this or re-use code: the code block is identical
     *        to that above, except for extraRateCoefficients instead of rateCoefficients.
     */
    if (extraRateCoefficients.size())
    {
        json_type& out = (*m_active)["rate_coefficients_extra"];
        out["labels"] = { "Ine.R.Coeff.", "Sup.R.Coeff.", "Description"};
        out["units"] = { "m^3/s", "m^3/s", ""};
        json_type& data = out["data"];
        for (const auto &rateCoeff : extraRateCoefficients)
        {
            std::stringstream ss;
            ss << *rateCoeff.collision;
            data.push_back( json_type{rateCoeff.inelastic, rateCoeff.superelastic, ss.str() } );
        }
    }
}

void JsonOutput::writeLookuptable(const Power &power, const SwarmParameters &swarmParameters) const
{
    /// \todo Handle isSimulationHF()
    json_type* out = m_root.contains("lookup_table")
                ? &m_root["lookup_table"]
                : nullptr;
    if (!out)
    {
        // make the section, and set up the labels and units
        out = &m_root["lookup_table"];
        /// \todo It would be nice to be able to add the quantities and units as pairs.
        (*out)["labels"] = { "RedField", "RedDif", "RedMob", "RedTow" "RedAtt","MeanE", "CharE", "EleTemp",
                      "DriftVelocity", "RelativePowerBalance" };
        (*out)["units"] = { "Td", "1/(m*s)", "1/(m*s*V)", "m^2" "m^2", "eV", "eV", "eV", "m/s", "1" };
    }
    assert(out);

    (*out)["data"].push_back( {
            m_workingConditions->reducedField(),
            swarmParameters.redDiffCoeff,
            swarmParameters.redMobCoeff,
            swarmParameters.redTownsendCoeff,
            swarmParameters.redAttCoeff,
            swarmParameters.meanEnergy,
            swarmParameters.characEnergy,
            swarmParameters.Te,
            swarmParameters.driftVelocity,
            power.relativeBalance * 100
    } );
}


} // namespace loki
