// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_TYPEDEFS_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_TYPEDEFS_H_

#include "Common/Constants.h"
#include "DynamicRupture/Misc.h"
#include "Kernels/Precision.h"
#include <Common/Executor.h>

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

template <Executor Executor>
struct FaultStresses;

template <Executor Executor>
struct TractionResults;

/**
 * Struct that contains all input stresses
 * normalStress in direction of the face normal, traction1, traction2 in the direction of the
 * respective tangential vectors
 */
template <>
struct FaultStresses<Executor::Host> {
  alignas(Alignment) real normalStress[ConvergenceOrder][misc::NumPaddedPoints] = {{}};
  alignas(Alignment) real traction1[ConvergenceOrder][misc::NumPaddedPoints] = {{}};
  alignas(Alignment) real traction2[ConvergenceOrder][misc::NumPaddedPoints] = {{}};
  alignas(Alignment) real fluidPressure[ConvergenceOrder][misc::NumPaddedPoints] = {{}};
};

/**
 * Struct that contains all traction results
 * traction1, traction2 in the direction of the respective tangential vectors
 */
template <>
struct TractionResults<Executor::Host> {
  alignas(Alignment) real traction1[ConvergenceOrder][misc::NumPaddedPoints] = {{}};
  alignas(Alignment) real traction2[ConvergenceOrder][misc::NumPaddedPoints] = {{}};
};

/**
 * Struct that contains all input stresses
 * normalStress in direction of the face normal, traction1, traction2 in the direction of the
 * respective tangential vectors
 */
template <>
struct FaultStresses<Executor::Device> {
  real normalStress[ConvergenceOrder] = {{}};
  real traction1[ConvergenceOrder] = {{}};
  real traction2[ConvergenceOrder] = {{}};
  real fluidPressure[ConvergenceOrder] = {{}};
};

/**
 * Struct that contains all traction results
 * traction1, traction2 in the direction of the respective tangential vectors
 */
template <>
struct TractionResults<Executor::Device> {
  real traction1[ConvergenceOrder] = {{}};
  real traction2[ConvergenceOrder] = {{}};
};

} // namespace seissol::dr

#endif // SEISSOL_SRC_DYNAMICRUPTURE_TYPEDEFS_H_
