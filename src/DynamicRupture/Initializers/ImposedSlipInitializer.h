#ifndef SEISSOL_IMPOSEDSLIPINITIALIZER_H
#define SEISSOL_IMPOSEDSLIPINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {
class ImposedSlipRatesInitializer;
}

/**
 * Derived initializer class for the ImposedSliprates friction law
 * Currently this is disabled, since Thomas is doing a lot of work on the master branch
 */
class seissol::dr::initializers::ImposedSlipRatesInitializer
    : public seissol::dr::initializers::BaseDRInitializer {
  public:
  using BaseDRInitializer::BaseDRInitializer;
};

#endif // SEISSOL_IMPOSEDSLIPINITIALIZER_H
