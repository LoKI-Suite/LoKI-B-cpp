//
// Created by daan on 2-5-19.
//

#include <Gas.h>

namespace loki {


    bool Gas::operator==(const Gas &other) {
        return name == other.name;
    }
}