#ifndef SEISSOL_IMPOSEDSLIPINITIALIZER_H
#define SEISSOL_IMPOSEDSLIPINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {

/**
 * Derived initializer class for the ImposedSliprates friction law
 * Currently this is disabled, since Thomas is doing a lot of work on the master branch
 */
class ImposedSlipRatesInitializer : public BaseDRInitializer {
  public:
  using BaseDRInitializer::BaseDRInitializer;
};
} // namespace seissol::dr::initializers

#endif // SEISSOL_IMPOSEDSLIPINITIALIZER_H
