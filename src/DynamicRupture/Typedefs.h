#ifndef DR_TYPEDEFS
#define DR_TYPEDEFS

#include "Common/Constants.h"
#include "DynamicRupture/Misc.h"
#include "Kernels/Precision.h"

namespace seissol::dr {

/**
 * Stores the P and S wave impedances for an element and its neighbor as well as the eta values from
 * Carsten Uphoff's dissertation equation (4.51)
 */
struct ImpedancesAndEta {
  real zp, zs, zpNeig, zsNeig, etaP, etaS, invEtaS, invZp, invZs, invZpNeig, invZsNeig;
};

/**
 * Stores the impedance matrices for an element and its neighbor for a poroelastic material.
 * This generalizes equation (4.51) from Carsten's thesis
 */
struct ImpedanceMatrices {
  alignas(Alignment) real impedance[tensor::Zplus::size()] = {};
  alignas(Alignment) real impedanceNeig[tensor::Zminus::size()] = {};
  alignas(Alignment) real eta[tensor::eta::size()] = {};
};

/**
 * Struct that contains all input stresses
 * normalStress in direction of the face normal, traction1, traction2 in the direction of the
 * respective tangential vectors
 */
struct FaultStresses {
  alignas(Alignment) real normalStress[ConvergenceOrder][misc::NumPaddedPoints] = {{}};
  alignas(Alignment) real traction1[ConvergenceOrder][misc::NumPaddedPoints] = {{}};
  alignas(Alignment) real traction2[ConvergenceOrder][misc::NumPaddedPoints] = {{}};
  alignas(Alignment) real fluidPressure[ConvergenceOrder][misc::NumPaddedPoints] = {{}};
};

/**
 * Struct that contains all traction results
 * traction1, traction2 in the direction of the respective tangential vectors
 */
struct TractionResults {
  alignas(Alignment) real traction1[ConvergenceOrder][misc::NumPaddedPoints] = {{}};
  alignas(Alignment) real traction2[ConvergenceOrder][misc::NumPaddedPoints] = {{}};
};

} // namespace seissol::dr

#endif
