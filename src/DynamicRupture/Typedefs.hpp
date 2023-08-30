#ifndef DR_TYPEDEFS
#define DR_TYPEDEFS

#include "DynamicRupture/Misc.h"
#include "Kernels/precision.hpp"

namespace seissol::dr {

/**
 * Stores the different types of friction laws
 * The values resemble the identifiers used in the old fortran implementation
 */
enum class FrictionLawType : unsigned int {
  NoFault = 0,
  LinearSlipWeakening = 16,
  LinearSlipWeakeningBimaterial = 6,
  RateAndStateAgingLaw = 3,
  RateAndStateSlipLaw = 4,
  RateAndStateFastVelocityWeakening = 103,
  ImposedSlipRatesYoffe = 33,
  ImposedSlipRatesGaussian = 34,
  RateAndStateVelocityWeakening = 7,
  RateAndStateAgingNucleation = 101,
  RateAndStateThermalProxy = 105,
};

/**
 * Stores the P and S wave impedances for an element and its neighbor as well as the eta values from
 * Carsten Uphoff's dissertation equation (4.51)
 */
struct ImpedancesAndEta {
  real zp, zs, zpNeig, zsNeig, etaP, etaS, invEtaS, invZp, invZs, invZpNeig, invZsNeig;
};

/**
 * Struct that contains all input stresses
 * normalStress in direction of the face normal, traction1, traction2 in the direction of the
 * respective tangential vectors
 */
struct FaultStresses {
  alignas(ALIGNMENT) real normalStress[CONVERGENCE_ORDER][misc::numPaddedPoints] = {{}};
  alignas(ALIGNMENT) real traction1[CONVERGENCE_ORDER][misc::numPaddedPoints] = {{}};
  alignas(ALIGNMENT) real traction2[CONVERGENCE_ORDER][misc::numPaddedPoints] = {{}};
};

/**
 * Struct that contains all traction results
 * traction1, traction2 in the direction of the respective tangential vectors
 */
struct TractionResults {
  alignas(ALIGNMENT) real traction1[CONVERGENCE_ORDER][misc::numPaddedPoints] = {{}};
  alignas(ALIGNMENT) real traction2[CONVERGENCE_ORDER][misc::numPaddedPoints] = {{}};
};

} // namespace seissol::dr

#endif
