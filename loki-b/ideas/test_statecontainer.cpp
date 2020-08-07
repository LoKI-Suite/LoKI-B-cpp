#include "ideas/StateContainer.h"
#include "ideas/PropertyFunctions.h"
#include <iostream>

int main()
{
    loki::experimental::GasMixture mixture;

    const auto& N2X0 = mixture.ensure_state("N2","0","X","0");
    mixture.ensure_state("N2","0","X","0","1");
    mixture.ensure_state("N2","0","X","0","2");
    mixture.ensure_state("N2","0","X","0","3");
    mixture.ensure_state("N2","0","X","0","4");
    mixture.ensure_state("N2","0","X","0","5");
    mixture.ensure_state("N2","0","X","0","6");
    mixture.ensure_state("N2","0","X","0","7");

    mixture.ensure_state("N2","0","X","1");
    mixture.ensure_state("N2","0","X","2");
    mixture.ensure_state("N2","0","X","3");

    loki::experimental::PropertyFunctions::rotationalDegeneracy(N2X0->children(),{});
    /// We need Gasmixture::getGas(string), that throws if the gas does not exists.
    loki::experimental::Gas& N2 = mixture.ensure_gas("N2");
    N2.rotationalConstant = 2.477501826053967e-4;
    N2.harmonicFrequency = 4.442724077107641875e+14;
    loki::experimental::PropertyFunctions::rigidRotorEnergy(N2X0->children(),{});
    loki::experimental::PropertyFunctions::boltzmannPopulation(N2X0->children(),{300.0});

    auto& N2X = *mixture.ensure_state("N2","0","X");
    loki::experimental::PropertyFunctions::harmonicOscillatorEnergy(N2X.children(),{});
    loki::experimental::PropertyFunctions::constantValue(N2X.children(), 1.0, loki::StatePropertyType::statisticalWeight);
    loki::experimental::PropertyFunctions::boltzmannPopulation(N2X.children(),{300.0});
    loki::experimental::PropertyFunctions::boltzmannPopulation(N2X.ensure_state("0")->children(),{300.0});

    mixture.ensure_state("O2","0","X","0","1");
    mixture.ensure_state("O2","0","X","0","3");
    mixture.ensure_state("O2","0","X","0","5");
    auto& O2_X_0_rot = *mixture.ensure_state("O2","0","X","0","7");

    N2.printStates(std::cout);
    return 0;
}
