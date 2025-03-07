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

template <typename T>
using PointData = T[misc::NumPaddedPoints];

/**
 * Stores the P and S wave impedances for an element and its neighbor as well as the eta values from
 * Carsten Uphoff's dissertation equation (4.51)
 */
template <bool Pointwise>
struct ImpedancesAndEta;

template <>
struct ImpedancesAndEta<false> {
  real zp_[1]{};
  real zs_[1]{};
  real zpNeig_[1]{};
  real zsNeig_[1]{};
  real etaP_[1]{};
  real etaS_[1]{};
  real invEtaS_[1]{};
  real invZp_[1]{};
  real invZs_[1]{};
  real invZpNeig_[1]{};
  real invZsNeig_[1]{};

#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real zp(int index) const { return zp_[0]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real zs(int index) const { return zs_[0]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real zpNeig(int index) const { return zpNeig_[0]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real zsNeig(int index) const { return zsNeig_[0]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real etaP(int index) const { return etaP_[0]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real etaS(int index) const { return etaS_[0]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real invEtaS(int index) const { return invEtaS_[0]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real invZp(int index) const { return invZp_[0]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real invZs(int index) const { return invZs_[0]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real invZpNeig(int index) const {
    return invZpNeig_[0];
  }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real invZsNeig(int index) const {
    return invZsNeig_[0];
  }
};

template <>
struct ImpedancesAndEta<true> {
  PointData<real> zp_{};
  PointData<real> zs_{};
  PointData<real> zpNeig_{};
  PointData<real> zsNeig_{};
  PointData<real> etaP_{};
  PointData<real> etaS_{};
  PointData<real> invEtaS_{};
  PointData<real> invZp_{};
  PointData<real> invZs_{};
  PointData<real> invZpNeig_{};
  PointData<real> invZsNeig_{};

#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real zp(int index) const { return zp_[index]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real zs(int index) const { return zs_[index]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real zpNeig(int index) const { return zpNeig_[index]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real zsNeig(int index) const { return zsNeig_[index]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real etaP(int index) const { return etaP_[index]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real etaS(int index) const { return etaS_[index]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real invEtaS(int index) const {
    return invEtaS_[index];
  }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real invZp(int index) const { return invZp_[index]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real invZs(int index) const { return invZs_[index]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real invZpNeig(int index) const {
    return invZpNeig_[index];
  }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real invZsNeig(int index) const {
    return invZsNeig_[index];
  }
};

/**
 * Stores the impedance matrices for an element and its neighbor for a poroelastic material.
 * This generalizes equation (4.51) from Carsten's thesis
 */
template <bool Pointwise>
struct ImpedanceMatrices;

template <>
struct ImpedanceMatrices<false> {
  alignas(Alignment) real impedance_[tensor::Zplus::size()] = {};
  alignas(Alignment) real impedanceNeig_[tensor::Zminus::size()] = {};
  alignas(Alignment) real eta_[tensor::eta::size()] = {};
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr const real* __restrict impedance(int index) const {
    return impedance_;
  }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr const real* __restrict impedanceNeig(int index) const {
    return impedanceNeig_;
  }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr const real* __restrict eta(int index) const {
    return eta_;
  }
};

template <>
struct ImpedanceMatrices<true> {
  alignas(Alignment) PointData<real[tensor::Zplus::size()]> impedance_ = {};
  alignas(Alignment) PointData<real[tensor::Zminus::size()]> impedanceNeig_ = {};
  alignas(Alignment) PointData<real[tensor::eta::size()]> eta_ = {};
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr const real* __restrict impedance(int index) const {
    return impedance_[index];
  }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr const real* __restrict impedanceNeig(int index) const {
    return impedanceNeig_[index];
  }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr const real* __restrict eta(int index) const {
    return eta_[index];
  }
};

template <bool Pointwise>
struct IsotropicWaveSpeeds;

template <>
struct IsotropicWaveSpeeds<false> {
  double density_[1];
  double pWaveVelocity_[1];
  double sWaveVelocity_[1];

#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr double density(int index) const { return density_[0]; }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr double pWaveVelocity(int index) const {
    return pWaveVelocity_[0];
  }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr double sWaveVelocity(int index) const {
    return sWaveVelocity_[0];
  }
};

template <>
struct IsotropicWaveSpeeds<true> {
  PointData<double> density_;
  PointData<double> pWaveVelocity_;
  PointData<double> sWaveVelocity_;

#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr double density(int index) const {
    return density_[index];
  }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr double pWaveVelocity(int index) const {
    return pWaveVelocity_[index];
  }
#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr double sWaveVelocity(int index) const {
    return sWaveVelocity_[index];
  }
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
  alignas(Alignment) PointData<real> normalStress[ConvergenceOrder] = {{}};
  alignas(Alignment) PointData<real> traction1[ConvergenceOrder] = {{}};
  alignas(Alignment) PointData<real> traction2[ConvergenceOrder] = {{}};
  alignas(Alignment) PointData<real> fluidPressure[ConvergenceOrder] = {{}};
};

/**
 * Struct that contains all traction results
 * traction1, traction2 in the direction of the respective tangential vectors
 */
template <>
struct TractionResults<Executor::Host> {
  alignas(Alignment) PointData<real> traction1[ConvergenceOrder] = {{}};
  alignas(Alignment) PointData<real> traction2[ConvergenceOrder] = {{}};
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
