#include "ideas/StateContainer.h"
#include <iostream>

int main()
{
    loki::Gas gas{"N2"};
    gas.add_state("0");
    gas.add_state("1","A");
    gas.add_state("1","B","1");
    gas.add_state("1","B","2","1");
    gas.add_state("1","B")->add_child("2","2");

    gas.printStates(std::cout);
    return 0;
}
