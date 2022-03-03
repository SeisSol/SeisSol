#ifndef DR_TYPEDEFS
#define DR_TYPEDEFS

#include "Kernels/precision.hpp"
#include "DynamicRupture/Misc.h"

namespace seissol::dr {

/**
 * Stores the different types of friction laws
 * The values resemble the identifiers used in the old fortran implementation
 */
enum class FrictionLawType : unsigned int {
  NoFault = 0,
  LinearSlipWeakening = 2,
  LinearSlipWeakeningBimaterial = 6,
  LinearSlipWeakeningForcedRuptureTime = 16,
  RateAndStateAgingLaw = 3,
  RateAndStateSlipLaw = 4,
  RateAndStateFastVelocityWeakening = 103,
  ImposedSlipRates = 33,
  RateAndStateVelocityWeakening = 7,
  RateAndStateAgingNucleation = 101,
};

/**
 * Stores the P and S wave impedances for an element and its neighbor as well as the eta values from
 * Uphoff equation (4.51)
 */
struct ImpedancesAndEta {
  real zp, zs, zpNeig, zsNeig, etaP, etaS, invEtaS, invZp, invZs, invZpNeig, invZsNeig;
};

/**
 * Struct that contains all input stresses
 */
struct FaultStresses {
  real normalStress[CONVERGENCE_ORDER][misc::numPaddedPoints] = {{}};
  real xyStress[CONVERGENCE_ORDER][misc::numPaddedPoints] = {{}};
  real xzStress[CONVERGENCE_ORDER][misc::numPaddedPoints] = {{}};
};

/**
 * Struct that contains all traction results
 */
struct TractionResults {
  real xyTraction[CONVERGENCE_ORDER][misc::numPaddedPoints] = {{}};
  real xzTraction[CONVERGENCE_ORDER][misc::numPaddedPoints] = {{}};
};

} // namespace seissol::dr

#endif
