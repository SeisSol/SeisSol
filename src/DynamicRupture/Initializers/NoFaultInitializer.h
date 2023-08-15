#ifndef SEISSOL_NOFAULTINITIALIZER_H
#define SEISSOL_NOFAULTINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {

/**
 * Derived initializer class for the NoFault friction law
 * does nothing in particular
 */
template <typename Config>
class NoFaultInitializer : public BaseDRInitializer<Config> {
  public:
  using BaseDRInitializer<Config>::BaseDRInitializer;
  using BaseDRInitializer<Config>::initializeFault;
};
} // namespace seissol::dr::initializers

#endif // SEISSOL_NOFAULTINITIALIZER_H
