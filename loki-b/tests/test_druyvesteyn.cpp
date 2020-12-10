#include <iostream>
#include <cmath>

#include "LoKI-B/Setup.h"
#include "LoKI-B/Simulation.h"

#define TOLERANCE 2e-5

auto json = R"json(
{
    "workingConditions": {
        "reducedField": 10.0,
        "electronTemperature": "0.03",
        "excitationFrequency": 0,
        "gasPressure": 133.32,
        "gasTemperature": 0,
        "electronDensity": 1e19,
        "chamberLength": 1.0,
        "chamberRadius": 1.0
    },
    "electronKinetics":    {
        "isOn": true,
        "eedfType": "boltzmann",
        "ionizationOperatorType": "usingSDCS",
        "growthModelType": "temporal",
        "includeEECollisions": false,
        "mixture": {
            "states": {
                "e": {
                    "particle": "e",
                    "charge": -1
                },
                "X(0)": {
                    "particle": "X",
                    "charge": 0,
                    "electronic": [
                        {
                            "e": "0",
                            "summary": "0"
                        }
                    ]
                },
                "X(*)": {
                    "particle": "X",
                    "charge": 0,
                    "electronic": [
                        {
                            "e": "*",
                            "summary": "*"
                        }
                    ]
                },
                "X(+)": {
                    "particle": "X",
                    "charge": 1,
                    "electronic": [
                        {
                            "e": "0",
                            "summary": "0"
                        }
                    ]
                }
            },
            "sets": [
                        {
                "processes": [
                    {
                        "reaction": {
                            "lhs": [
                                { "state": "e", "count": 1 },
                                { "state": "X(0)", "count": 1 }
                        ],
                        "rhs": [
                            { "state": "e", "count": 1 },
                            { "state": "X(0)", "count": 1 }
                        ],
                        "reversible": false,
                        "type_tags": ["Elastic"]
                        },
                        "parameters": ["m/M = 0.00025"],
                        "reference": [
                            "Analytical test, constant cross section"
                        ],
                        "labels": ["Energy", "Cross section"],
                        "units": ["eV", "m2"],
                        "threshold": 0,
                        "data": [
                            [0.000000e+0, 1.000000e-19],
                            [1.000000e+3, 1.000000e-19]
                        ]
                    }
                ]
            }
            ]
        },
        "gasProperties": {
            "_mass": "test_masses.txt",
            "mass": [
                { "name": "X", "mass": 1.66e-27 }
            ],
            "fraction": [
                "X = 1"
            ]
        },
        "stateProperties": {
            "energy": [
                "X(0) = 0.0"
            ],
            "population": [
                "X(0) = 1.0"
            ],
            "statisticalWeight": [
                "X(0) = 1"
            ]
        },
        "numerics": {
            "energyGrid": {
                "maxEnergy": 1,
                "cellNumber": 1000,
                "smartGrid": {
                    "minEedfDecay": 19,
                    "maxEedfDecay": 21,
                    "updateFactor": 0.05
                }
            },
            "maxPowerBalanceRelError": 1e-9,
            "nonLinearRoutines": {
                "algorithm": "mixingDirectSolutions",
                "mixingParameter": 0.7,
                "maxEedfRelError": 1e-9
            }
        }
    },
    "gui": {
        "isOn": false,
        "refreshFrequency": 1
    },
    "output": {
        "isOn": true,
        "folder": "simulation",
        "dataFiles": [
            "eedf",
            "swarmParameters",
            "rateCoefficients",
            "powerBalance",
            "lookUpTable"
        ]
    }
}
)json"_json;

double druyvesteyn(const double eps, const double e)
{
    const double pi = 3.141592653589793238463;
    const double Gamma_1_4 = 3.625609908221908311931;
    const double Gamma_5_4 = Gamma_1_4 / 4;
    const double Gamma_3_4 = pi * std::sqrt(2) / Gamma_1_4;
    const double b_1 = std::pow(Gamma_5_4, 1.5) / std::pow(Gamma_3_4, 2.5);
    const double b_2 = Gamma_5_4 / Gamma_3_4;
    return 2 * b_1 * std::exp(-std::pow(b_2*e/eps, 2)) / std::pow(eps, 1.5);
}

void checkRMSE(const loki::Grid &grid,
    const loki::Vector &eedf,
    const loki::WorkingConditions &wc,
    const loki::Power &power,
    const std::vector<loki::EedfGas *> &gases,
    const loki::SwarmParameters &swarmParameters,
    const std::vector<loki::RateCoefficient> &rateCoefficients,
    const std::vector<loki::RateCoefficient> &extraRateCoefficients,
    const loki::Vector &firstAnisotropy)
{
    using namespace loki;
    if (grid.getCells().size() != eedf.size())
    {
        throw std::runtime_error("error");
    }

    const double eps = swarmParameters.meanEnergy;
    double se = 0.0;
    for (unsigned i = 0; i < grid.getCells().size(); i++)
    {
        const double e = grid.getCells()[i];
        const double eedf_d = druyvesteyn(eps, e);
        se += std::pow(eedf[i] - eedf_d, 2);
        std::cout << e << "\t" << eedf[i] << "\t" << eedf_d << std::endl;
    }
    const double rmse = std::sqrt(se / grid.getCells().size());
    if (rmse > TOLERANCE)
    {
        throw std::runtime_error(std::string("RMSE = ") + std::to_string(rmse) +
            " exceeds tolerance");
    }
}

int main(int argc, char **argv)
{
    try
    {
        std::unique_ptr<loki::Simulation> simulation(new loki::Simulation(json));
        simulation->m_obtainedResults.addListener(checkRMSE);
        simulation->run();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}
