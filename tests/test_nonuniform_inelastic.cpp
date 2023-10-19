/** \file
 *
 *  Unit tests for nonuniform inelastic operator
 *
 *  \author Jop Hendrikx
 *  \date   October 2023
 */

#include "LoKI-B/Grid.h"
#include "LoKI-B/EedfUtilities.h"
#include "LoKI-B/EedfCollisions.h"
#include "tests/TestUtilities.h"
#include "source/Operators.cpp"
#include <iostream>
#include <fstream>

auto json = R"json(
{
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
                "type_tags": ["Excitation"]
                },
                "parameters": ["m/M = 0.000025"],
                "reference": [
                    "Analytical test, constant cross section"
                ],
                "labels": ["Energy", "Cross section"],
                "units": ["eV", "m2"],
                "threshold": 9.99,
                "data": [
                    [9.990000e+0, 0],
                    [10.0, 1e6],
                    [10.010000e+0, 0]
                ]
            },
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
                    [-0.1,0],
                    [4.99,0], 
                    [5.0, 0.2],
                    [10.0, 0.2]
                ]
            },
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
                "threshold": 0.98,
                "data": [
                    [0.98,0],
                    [0.9990000e+0, 1e25],
                    [1.100000e+0, 1e25],
                    [1.11, 0]
                ]
            }  
        ]
    },
    "gasProperties": {
        "_mass": "test_masses.txt",
        "mass": [
            { "name": "X", "mass": 1 }
        ],
        "fraction": [
            "X = 1"
        ],
        "_harmonicFrequency": "test_harmonicFrequency.txt",
        "harmonicFrequency": [
            { "name": "X", "harmonicFrequency": 1 }
        ],
        "_anharmonicFrequency": "test_anharmonicFrequency.txt",
        "anharmonicFrequency": [
            { "name": "X", "anharmonicFrequency": 1 }
        ],
        "_electricQuadrupoleMoment": "test_electricQuadrupoleMoment.txt",
        "electricQuadrupoleMoment": [
            { "name": "X", "electricQuadrupoleMoment": 1 }
        ],
        "_rotationalConstant": "test_RotationalConstant.txt",
        "rotationalConstant": [
            { "name": "X", "rotationalConstant": 1 }
        ],
        "_OPBParameter": "test_OPBParameter.txt",
        "OPBParameter": [
            { "name": "X", "OPBParameter": 1 }
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
    }
    
}
)json"_json;

void write_eedf(std::ostream& os, const loki::Grid& grid, const loki::Vector& eedf)
{
    for (loki::Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        os << grid.getCells()(k) << '\t' << eedf(k) << std::endl;
    }
}

int main()
{
    using namespace loki;
    
    int nCells = 3000;
    double uMax = 30;
    double Tg = 00;
    const double eon = 10;
    const double won = 0;

    Grid::Vector nodeDistribution = Vector::LinSpaced(nCells + 1, 0.0, 1.0);
    Grid grid1(nodeDistribution, uMax, false);
    Grid grid2(nCells, uMax);

    ElectronKineticsSetup kineticSetup;
    WorkingConditionsSetup conditionsSetup;
    conditionsSetup.gasTemperature = Tg;
    WorkingConditions workConditions(conditionsSetup);
    
    EedfMixture mixture1(&grid1, json, &workConditions);
    EedfMixture mixture2(&grid2, json, &workConditions);
    
    InelasticOperator inelasticOperator1(grid1);
    InelasticOperator inelasticOperator2(grid2);
   
    SparseMatrix M1(nCells, nCells);
    Matrix boltzmann1 = M1.toDense(); 

    SparseMatrix M2(nCells, nCells);
    Matrix boltzmann2 = M2.toDense(); 

    
    inelasticOperator2.evaluateInelasticOperators(grid2,mixture2);
    inelasticOperator1.evaluateInelasticOperators(grid1,mixture1);
    test_expr(inelasticOperator1.inelasticMatrix.isApprox(inelasticOperator2.inelasticMatrix));


    // Elastic part
    Vector elasticCrossSection = Vector::Ones(nCells+1);
    ElasticOperator elasticOperator;
    SparseMatrix Melastic1(nCells, nCells);
    SparseMatrix Melastic2(nCells, nCells);
    
    elasticOperator.evaluate(grid1, elasticCrossSection, Tg, Melastic1);
    elasticOperator.evaluate(grid2, elasticCrossSection, Tg, Melastic2);
    test_expr(Melastic1.isApprox(Melastic2));
    
    // Field part
    Vector fieldCrossSection = Vector::Ones(nCells+1);
    FieldOperator fieldOperator1(grid1);
    FieldOperator fieldOperator2(grid2);
    SparseMatrix Mfield1(nCells, nCells);
    SparseMatrix Mfield2(nCells, nCells);
    fieldOperator1.evaluate(grid1, fieldCrossSection, eon, won, Mfield1);
    fieldOperator2.evaluate(grid2, fieldCrossSection, eon, won, Mfield2);

    boltzmann1 = inelasticOperator1.inelasticMatrix + Melastic1 + Mfield1;
    Vector eedf1 = Vector::Zero(nCells);
    eedf1[0] = 1.;

    boltzmann1.row(0) = grid1.getCells().cwiseSqrt().cwiseProduct(grid1.duCells());
    LinAlg::hessenberg(boltzmann1.data(), eedf1.data(), grid1.nCells());


    Vector f = Vector::Zero(nCells);
    Vector fIhg = Vector::Zero(nCells);

    double refU = 0;
    double sigma = 1e-1;
    double Q0 = 1e-0;
    // double mM = 2.5e-5;
    // double e = Constant::electronCharge;
    double u0 = 10;
    // double ureac = 10;
    // double W = e/ std::sqrt(6*mM) * eon/Q0;

    Grid::Index findIndex = (std::upper_bound(grid1.getCells().begin(),grid1.getCells().end(), u0) - grid1.getCells().begin());
    
    for (Grid::Index k = 0; k < grid1.nCells(); ++k)
        {
            // Analytical solution
            f[k] = std::exp(-3*Q0*(std::pow(grid1.getCell(k),2) - std::pow(refU,2)) /std::pow(eon,2));
            fIhg[k] = 6*sigma / (std::pow(eon,2) * Q0) * std::exp(-3*std::pow(Q0 * refU /eon,2))*(+std::expint(3*std::pow(Q0*grid1.getCell(k)/eon,2)) - std::expint(3*std::pow(Q0*u0/eon,2)));

            // Solution from Peter Koelman 
            // f[k] = std::exp(-1*std::pow(u0/W,2) / 2 * (std::pow(grid1.getCell(k)/u0,2) - std::pow(refU/u0,2)));
            // fihg[k] = sigma / 4 / mM / Q0 * std::pow(ureac/W,2) * std::exp(-std::pow(ureac/W,2)/2)*(std::expint(std::pow(grid1.getCell(k)/W,2)/2) - std::expint(std::pow(ureac/W,2)/2));
            if (k > findIndex)
            {
                f[k] *= (1+fIhg[k]);

            } 
            
        }
    f /= f.dot(grid1.getCells().cwiseSqrt().cwiseProduct(grid1.duCells())); 
    
    double relativeError = ((eedf1-f).cwiseQuotient(f).cwiseAbs()).mean();
    std::cout << relativeError << std::endl;

    std::ofstream ofs("eedfinelastic.dat");
    write_eedf(ofs,grid1,eedf1);

    test_report;
    return nerrors;
}
