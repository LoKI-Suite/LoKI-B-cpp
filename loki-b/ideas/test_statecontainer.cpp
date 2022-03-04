#include "ideas/StateContainer.h"
#include "ideas/PropertyFunctions.h"
#include <iostream>

int main()
{
    using namespace loki::experimental;

    GasMixture mixture;

    StateRoot* N2_root = mixture.ensure_state("N2");
    StateCharge* N20 = N2_root->ensure_state("0");
    StateElectronic* N20X = N20->ensure_state("X");
    StateVibrational* N20X0 = N20X->ensure_state("0");
    /* StateRotational* N20X01 =*/ N20X0->ensure_state("1");
    // and J=2...6
    N20X0->ensure_state("2");
    N20X0->ensure_state("3");
    N20X0->ensure_state("4");
    N20X0->ensure_state("5");
    N20X0->ensure_state("6");
    // and J=7, in one go, starting from the mixture: arguments are gas, charge, elec, vib, rot
    mixture.ensure_state("N2","0","X","0","7");
    // three more vibrational states:
    mixture.ensure_state("N2","0","X","1");
    mixture.ensure_state("N2","0","X","2");
    mixture.ensure_state("N2","0","X","3");
    // create N2(+)
    /* StateCharge* N21 =*/ N2_root->ensure_state("1");

    PropertyFunctions::rotationalDegeneracy(N20X0->children(),{});
    /// \todo We need Gasmixture::getGas(string), that throws if the gas does not exists.
    Gas& N2 = mixture.ensure_gas("N2");
    N2.rotationalConstant = 2.477501826053967e-4;
    N2.harmonicFrequency = 4.442724077107641875e+14;
    PropertyFunctions::rigidRotorEnergy(N20X0->children(),{});
    PropertyFunctions::boltzmannPopulation(N20X0->children(),{300.0});

    N20->population() = 1.0;
    N20X->population() = 1.0;
    PropertyFunctions::harmonicOscillatorEnergy(N20X->children(),{});
    PropertyFunctions::constantValue(N20X->children(), 1.0, loki::StatePropertyType::statisticalWeight);
    PropertyFunctions::boltzmannPopulation(N20X->children(),{300.0});
    PropertyFunctions::boltzmannPopulation(N20X->ensure_state("0")->children(),{300.0});

    mixture.ensure_state("O2","0","X","0","1");
    mixture.ensure_state("O2","0","X","0","3");
    mixture.ensure_state("O2","0","X","0","5");
    mixture.ensure_state("O2")->population() = 1.0;
    mixture.ensure_state("O2","0")->population() = 1.0;
    mixture.ensure_state("O2","0","X")->population() = 1.0;
    mixture.ensure_state("O2","0","X","0")->population() = 1.0;
    mixture.ensure_state("O2","0","X","0","1")->population() = 1.0;
    mixture.ensure_state("O2","0","X","0","3")->population() = 0.0;
    mixture.ensure_state("O2","0","X","0","5")->population() = 0.0;

    mixture.print(std::cout);

    mixture.checkPopulations();
    mixture.evaluateStateDensities();

    apply(mixture.ensure_gas("O2").states(), [](const StateData& st) { std::cout << st.density() << std::endl; }, true );
    apply(mixture.ensure_gas("O2").states(), [](const auto& st) { std::cout << std::string(st.level(),' ') << st.str() << std::endl; }, true );
    return 0;
}
