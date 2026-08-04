// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_TYPEDEFS_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_TYPEDEFS_H_

#include "Alignment.h"
#include "Common/Constants.h"
#include "Common/Executor.h"
#include "DynamicRupture/Misc.h"
#include "Kernels/Precision.h"

namespace seissol::dr {

/**
 * Stores the P and S wave impedances for an element and its neighbor as well as the eta values from
 * Carsten Uphoff's dissertation equation (4.51)
 */
struct ImpedancesAndEta {
  real zp{};
  real zs{};
  real zpNeig{};
  real zsNeig{};
  real etaP{};
  real etaS{};
  real invEtaS{};
  real invZp{};
  real invZs{};
  real invZpNeig{};
  real invZsNeig{};
};

/**
 * Stores the impedance matrices for an element and its neighbor for a poroelastic material.
 * This generalizes equation (4.51) from Carsten's thesis
 */
struct ImpedanceMatrices {
  alignas(Alignment) real impedance[tensor::Zplus::size()] = {};
  alignas(Alignment) real impedanceNeig[tensor::Zminus::size()] = {};
  alignas(Alignment) real eta[tensor::eta::size()] = {};
  /**
   * Maps a fault-local traction difference to the difference of the stress components which do not
   * take part in the fault-normal Riemann problem:
   *
   *   [d sigma_ss; d sigma_dd; d sigma_sd]
   *       = lateralStress * [d sigma_nn; d sigma_ns; d sigma_nd]
   *
   * For anisotropy this is C[{ss,dd,sd},{nn,ns,nd}] * Gamma^-1 with the Christoffel matrix Gamma of
   * the fault normal direction; for poroelasticity the traction carries a fourth component (the
   * fluid pressure) and the matrix is 3 x 4. Dense and column major, lateralStress[col * 3 + row].
   *
   * Only needed by the fault receiver output, and only filled for the materials that need the
   * matrix form -- for an isotropic elastic material the single relevant entry is
   * lambda / (lambda + 2 mu) = 1 - 2 (cs/cp)^2, which the output computes from the wave speeds.
   */
  alignas(Alignment) real lateralStress[3 * tensor::Zminus::Shape[0]] = {};
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
  alignas(Alignment) real normalStress[misc::TimeSteps][misc::NumPaddedPoints] = {{}};
  alignas(Alignment) real traction1[misc::TimeSteps][misc::NumPaddedPoints] = {{}};
  alignas(Alignment) real traction2[misc::TimeSteps][misc::NumPaddedPoints] = {{}};
  alignas(Alignment) real fluidPressure[misc::TimeSteps][misc::NumPaddedPoints] = {{}};
};

/**
 * Struct that contains all traction results
 * normalStress in direction of the face normal, traction1, traction2 in the direction of the
 * respective tangential vectors.
 *
 * normalStress is the fault-normal traction *after* the friction solve. It differs from
 * FaultStresses::normalStress only if the impedance couples the fault-normal and the tangential
 * directions, which is the case for anisotropic materials. It is therefore seeded with the trial
 * value by common::initializeTractionResults, and a friction law only has to overwrite it if it
 * actually changes the normal stress.
 */
template <>
struct TractionResults<Executor::Host> {
  alignas(Alignment) real normalStress[misc::TimeSteps][misc::NumPaddedPoints] = {{}};
  alignas(Alignment) real traction1[misc::TimeSteps][misc::NumPaddedPoints] = {{}};
  alignas(Alignment) real traction2[misc::TimeSteps][misc::NumPaddedPoints] = {{}};
};

/**
 * Struct that contains all input stresses
 * normalStress in direction of the face normal, traction1, traction2 in the direction of the
 * respective tangential vectors
 */
template <>
struct FaultStresses<Executor::Device> {
  real normalStress[misc::TimeSteps] = {{}};
  real traction1[misc::TimeSteps] = {{}};
  real traction2[misc::TimeSteps] = {{}};
  real fluidPressure[misc::TimeSteps] = {{}};
};

/**
 * Struct that contains all traction results
 * normalStress in direction of the face normal, traction1, traction2 in the direction of the
 * respective tangential vectors. See TractionResults<Executor::Host> for the semantics of
 * normalStress.
 */
template <>
struct TractionResults<Executor::Device> {
  real normalStress[misc::TimeSteps] = {{}};
  real traction1[misc::TimeSteps] = {{}};
  real traction2[misc::TimeSteps] = {{}};
};

} // namespace seissol::dr

#endif // SEISSOL_SRC_DYNAMICRUPTURE_TYPEDEFS_H_
