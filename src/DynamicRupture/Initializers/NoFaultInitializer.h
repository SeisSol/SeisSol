#ifndef SEISSOL_NOFAULTINITIALIZER_H
#define SEISSOL_NOFAULTINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {
class NoFaultInitializer; // NoFaultFL0
}

// does not need implementation, inherits everything from BaseDR
class seissol::dr::initializers::NoFaultInitializer
    : public seissol::dr::initializers::BaseDRInitializer {
  public:
  using BaseDRInitializer::BaseDRInitializer;
};

#endif // SEISSOL_NOFAULTINITIALIZER_H
