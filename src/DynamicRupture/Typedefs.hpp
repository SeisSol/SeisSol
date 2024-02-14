#ifndef DR_TYPEDEFS
#define DR_TYPEDEFS

#include "DynamicRupture/Misc.h"
#include "Kernels/precision.hpp"

namespace seissol::dr {

/**
 * Stores the P and S wave impedances for an element and its neighbor as well as the eta values from
 * Carsten Uphoff's dissertation equation (4.51)
 */
struct ImpedancesAndEta {
  real zp, zs, zpNeig, zsNeig, etaP, etaS, invEtaS, invZp, invZs, invZpNeig, invZsNeig;
#if defined USE_DAMAGEDELASTIC
  real lambda0P, mu0P, gammaRP, xi0P, rho0P;
  real lambda0M, mu0M, gammaRM, xi0M, rho0M;
  real csOcpTZsOZp, csOcpTZsOZpNeig;
  VrtxCoords faultN, faultT1, faultT2;
#endif
};

/**
 * Stores the impedance matrices for an element and its neighbor for a poroelastic material.
 * This generalizes equation (4.51) from Carsten's thesis
 */
struct ImpedanceMatrices {
  alignas(ALIGNMENT) real impedance[tensor::Zplus::size()] = {};
  alignas(ALIGNMENT) real impedanceNeig[tensor::Zminus::size()] = {};
  alignas(ALIGNMENT) real eta[tensor::eta::size()] = {};
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
  alignas(ALIGNMENT) real fluidPressure[CONVERGENCE_ORDER][misc::numPaddedPoints] = {{}};
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
