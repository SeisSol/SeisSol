#ifndef SEISSOL_NOFAULTINITIALIZER_H
#define SEISSOL_NOFAULTINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {

/**
 * Derived initializer class for the NoFault friction law
 * does nothing in particular
 */
class NoFaultInitializer : public BaseDRInitializer {
  public:
  using BaseDRInitializer::BaseDRInitializer;
  using BaseDRInitializer::initializeFault;
};
} // namespace seissol::dr::initializers

#endif // SEISSOL_NOFAULTINITIALIZER_H
