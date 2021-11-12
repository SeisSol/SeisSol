#ifndef SEISSOL_NOFAULTINITIALIZER_H
#define SEISSOL_NOFAULTINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {
class NoFaultInitializer; // NoFaultFL0
}

/**
 * Derived initializer class for the NoFault friction law
 * does nothing in particular
 */
class seissol::dr::initializers::NoFaultInitializer
    : public seissol::dr::initializers::BaseDRInitializer {
  public:
  using BaseDRInitializer::BaseDRInitializer;
  using BaseDRInitializer::initializeFault;
};

#endif // SEISSOL_NOFAULTINITIALIZER_H
