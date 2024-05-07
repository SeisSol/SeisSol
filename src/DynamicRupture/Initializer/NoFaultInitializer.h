#ifndef SEISSOL_NOFAULTINITIALIZER_H
#define SEISSOL_NOFAULTINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializer {

/**
 * Derived initializer class for the NoFault friction law
 * does nothing in particular
 */
class NoFaultInitializer : public BaseDRInitializer {
  public:
  using BaseDRInitializer::BaseDRInitializer;
  using BaseDRInitializer::initializeFault;
};
} // namespace seissol::dr::initializer

#endif // SEISSOL_NOFAULTINITIALIZER_H
