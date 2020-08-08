#include "ideas/StateContainer.h"
#include "ideas/PropertyFunctions.h"
#include <iostream>

int main()
{
    loki::experimental::GasMixture mixture;

    const auto& N20X0 = *mixture.ensure_state("N2","0","X","0");
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

    loki::experimental::PropertyFunctions::rotationalDegeneracy(N20X0.children(),{});
    /// We need Gasmixture::getGas(string), that throws if the gas does not exists.
    loki::experimental::Gas& N2 = mixture.ensure_gas("N2");
    N2.rotationalConstant = 2.477501826053967e-4;
    N2.harmonicFrequency = 4.442724077107641875e+14;
    loki::experimental::PropertyFunctions::rigidRotorEnergy(N20X0.children(),{});
    loki::experimental::PropertyFunctions::boltzmannPopulation(N20X0.children(),{300.0});

    auto& N20 = *mixture.ensure_state("N2","0");
    N20.info().population = 1.0;
    auto& N20X = *mixture.ensure_state("N2","0","X");
    N20X.info().population = 1.0;
    loki::experimental::PropertyFunctions::harmonicOscillatorEnergy(N20X.children(),{});
    loki::experimental::PropertyFunctions::constantValue(N20X.children(), 1.0, loki::StatePropertyType::statisticalWeight);
    loki::experimental::PropertyFunctions::boltzmannPopulation(N20X.children(),{300.0});
    loki::experimental::PropertyFunctions::boltzmannPopulation(N20X.ensure_state("0")->children(),{300.0});

    mixture.ensure_state("O2","0","X","0","1");
    mixture.ensure_state("O2","0","X","0","3");
    mixture.ensure_state("O2","0","X","0","5");
    mixture.ensure_state("O2","0")->info().population = 1.0;
    mixture.ensure_state("O2","0","X")->info().population = 1.0;
    mixture.ensure_state("O2","0","X","0")->info().population = 1.0;
    mixture.ensure_state("O2","0","X","0","1")->info().population = 1.0;
    mixture.ensure_state("O2","0","X","0","3")->info().population = 0.0;
    mixture.ensure_state("O2","0","X","0","5")->info().population = 0.0;

    mixture.print(std::cout);

    mixture.checkPopulations();
    return 0;
}
