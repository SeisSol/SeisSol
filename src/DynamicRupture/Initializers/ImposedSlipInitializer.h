#ifndef SEISSOL_IMPOSEDSLIPINITIALIZER_H
#define SEISSOL_IMPOSEDSLIPINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {
class ImposedSlipRatesInitializer; // imposed slip rates on boundary
}

/*
 * nucleationStressInFaultCS initialized which is used to impose slip rates on the fault surface
 */
// currently disabled, as Thomas is doing a lot of work on this FL on the master branch
class seissol::dr::initializers::ImposedSlipRatesInitializer
    : public seissol::dr::initializers::BaseDRInitializer {
  public:
  using BaseDRInitializer::BaseDRInitializer;
};

#endif // SEISSOL_IMPOSEDSLIPINITIALIZER_H
