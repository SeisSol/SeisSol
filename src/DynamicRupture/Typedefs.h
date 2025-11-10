// SPDX-FileCopyrightText: 2021 SeisSol Group
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
#include <Alignment.h>
#include <Common/Executor.h>
#include <GeneratedCode/tensor.h>

namespace seissol::dr {

/**
 * Stores the P and S wave impedances for an element and its neighbor as well as the eta values from
 * Carsten Uphoff's dissertation equation (4.51)
 */
template <typename RealT>
struct ImpedancesAndEta {
  RealT zp, zs, zpNeig, zsNeig, etaP, etaS, invEtaS, invZp, invZs, invZpNeig, invZsNeig;
};

/**
 * Stores the impedance matrices for an element and its neighbor for a poroelastic material.
 * This generalizes equation (4.51) from Carsten's thesis
 */
template <typename Cfg>
struct ImpedanceMatrices {
  alignas(Alignment) Real<Cfg> impedance[tensor::Zplus<Cfg>::size()] = {};
  alignas(Alignment) Real<Cfg> impedanceNeig[tensor::Zminus<Cfg>::size()] = {};
  alignas(Alignment) Real<Cfg> eta[tensor::eta<Cfg>::size()] = {};
};

template <typename Cfg, Executor Executor>
struct FaultStresses;

template <typename Cfg, Executor Executor>
struct TractionResults;

/**
 * Struct that contains all input stresses
 * normalStress in direction of the face normal, traction1, traction2 in the direction of the
 * respective tangential vectors
 */
template <typename Cfg>
struct FaultStresses<Cfg, Executor::Host> {
  alignas(Alignment) Real<Cfg> normalStress[misc::TimeSteps<Cfg>][misc::NumPaddedPoints<Cfg>] = {
      {}};
  alignas(Alignment) Real<Cfg> traction1[misc::TimeSteps<Cfg>][misc::NumPaddedPoints<Cfg>] = {{}};
  alignas(Alignment) Real<Cfg> traction2[misc::TimeSteps<Cfg>][misc::NumPaddedPoints<Cfg>] = {{}};
  alignas(Alignment) Real<Cfg> fluidPressure[misc::TimeSteps<Cfg>][misc::NumPaddedPoints<Cfg>] = {
      {}};
};

/**
 * Struct that contains all traction results
 * traction1, traction2 in the direction of the respective tangential vectors
 */
template <typename Cfg>
struct TractionResults<Cfg, Executor::Host> {
  alignas(Alignment) Real<Cfg> traction1[misc::TimeSteps<Cfg>][misc::NumPaddedPoints<Cfg>] = {{}};
  alignas(Alignment) Real<Cfg> traction2[misc::TimeSteps<Cfg>][misc::NumPaddedPoints<Cfg>] = {{}};
};

/**
 * Struct that contains all input stresses
 * normalStress in direction of the face normal, traction1, traction2 in the direction of the
 * respective tangential vectors
 */
template <typename Cfg>
struct FaultStresses<Cfg, Executor::Device> {
  Real<Cfg> normalStress[misc::TimeSteps<Cfg>] = {{}};
  Real<Cfg> traction1[misc::TimeSteps<Cfg>] = {{}};
  Real<Cfg> traction2[misc::TimeSteps<Cfg>] = {{}};
  Real<Cfg> fluidPressure[misc::TimeSteps<Cfg>] = {{}};
};

/**
 * Struct that contains all traction results
 * traction1, traction2 in the direction of the respective tangential vectors
 */
template <typename Cfg>
struct TractionResults<Cfg, Executor::Device> {
  Real<Cfg> traction1[misc::TimeSteps<Cfg>] = {{}};
  Real<Cfg> traction2[misc::TimeSteps<Cfg>] = {{}};
};

} // namespace seissol::dr

#endif // SEISSOL_SRC_DYNAMICRUPTURE_TYPEDEFS_H_
