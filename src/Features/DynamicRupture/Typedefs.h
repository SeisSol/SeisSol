// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_TYPEDEFS_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_TYPEDEFS_H_

#include "Alignment.h"
#include "Core/Common/Constants.h"
#include "Core/Common/Executor.h"
#include "Core/Common/Real.h"
#include "Core/IO/Datatype/Datatype.h"
#include "Core/IO/Datatype/Inference.h"
#include "Features/DynamicRupture/Misc.h"
#include "GeneratedCode/tensor.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace seissol::dr {

/**
 * Stores the P and S wave impedances for an element and its neighbor as well as the eta values from
 * Carsten Uphoff's dissertation equation (4.51)
 */
template <typename Cfg>
struct ImpedancesAndEta {
  Real<Cfg> zp{};
  Real<Cfg> zs{};
  Real<Cfg> zpNeig{};
  Real<Cfg> zsNeig{};
  Real<Cfg> etaP{};
  Real<Cfg> etaS{};
  Real<Cfg> invEtaS{};
  Real<Cfg> invZp{};
  Real<Cfg> invZs{};
  Real<Cfg> invZpNeig{};
  Real<Cfg> invZsNeig{};
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
  alignas(Alignment) Real<Cfg> lateralStress[3 * tensor::Zminus<Cfg>::Shape[0]] = {};
};

template <typename Cfg, Executor Executor>
struct FaultStresses;

template <typename Cfg, Executor Executor>
struct TractionResults;

template <typename Cfg, Executor Executor>
struct ImposedState;

/**
 * Struct that contains all input stresses
 * normalStress in direction of the face normal, traction1, traction2 in the direction of the
 * respective tangential vectors
 */
template <typename Cfg>
struct FaultStresses<Cfg, Executor::Host> {
  using real = Real<Cfg>;

  alignas(Alignment) real normalStress[misc::NumPaddedPoints<Cfg>]{};
  alignas(Alignment) real traction1[misc::NumPaddedPoints<Cfg>]{};
  alignas(Alignment) real traction2[misc::NumPaddedPoints<Cfg>]{};
  alignas(Alignment) real fluidPressure[misc::NumPaddedPoints<Cfg>]{};
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
template <typename Cfg>
struct TractionResults<Cfg, Executor::Host> {
  using real = Real<Cfg>;

  alignas(Alignment) real normalStress[misc::NumPaddedPoints<Cfg>]{};
  alignas(Alignment) real traction1[misc::NumPaddedPoints<Cfg>]{};
  alignas(Alignment) real traction2[misc::NumPaddedPoints<Cfg>]{};
};

/**
 * Accumulator for the imposed state. Used after every internal timestep.
 */
template <typename Cfg>
struct ImposedState<Cfg, Executor::Host> {
  using real = Real<Cfg>;

  alignas(Alignment) real plus[misc::NumQuantities<Cfg>][misc::NumPaddedPoints<Cfg>]{};
  alignas(Alignment) real minus[misc::NumQuantities<Cfg>][misc::NumPaddedPoints<Cfg>]{};
};

/**
 * Struct that contains all input stresses
 * normalStress in direction of the face normal, traction1, traction2 in the direction of the
 * respective tangential vectors
 */
template <typename Cfg>
struct FaultStresses<Cfg, Executor::Device> {
  using real = Real<Cfg>;

  real normalStress{};
  real traction1{};
  real traction2{};
  real fluidPressure{};
};

/**
 * Struct that contains all traction results
 * normalStress in direction of the face normal, traction1, traction2 in the direction of the
 * respective tangential vectors. See TractionResults<Cfg, Executor::Host> for the semantics of
 * normalStress.
 */
template <typename Cfg>
struct TractionResults<Cfg, Executor::Device> {
  using real = Real<Cfg>;

  real normalStress{};
  real traction1{};
  real traction2{};
};

/**
 * Accumulator for the imposed state. Used after every internal timestep.
 */
template <typename Cfg>
struct ImposedState<Cfg, Executor::Device> {
  using real = Real<Cfg>;

  real plus[misc::NumQuantities<Cfg>]{};
  real minus[misc::NumQuantities<Cfg>]{};
};

} // namespace seissol::dr

namespace seissol {

struct DRFaceInformation {
  std::size_t meshFace{};
  std::size_t copy{};
  std::uint8_t plusSide{};
  std::uint8_t minusSide{};
  std::uint8_t faceRelation{};
  bool plusSideOnThisRank{};
};

template <typename Cfg>
struct DRGodunovData {
  Real<Cfg> dataTinvT[seissol::tensor::TinvT<Cfg>::size()]{};
  Real<Cfg> tractionPlusMatrix[seissol::tensor::tractionPlusMatrix<Cfg>::size()]{};
  Real<Cfg> tractionMinusMatrix[seissol::tensor::tractionMinusMatrix<Cfg>::size()]{};
  // When integrating quantities over the fault
  // we need to integrate over each physical element.
  // The integration is effectively done in the reference element, and the scaling factor of
  // the transformation, the surface Jacobian (e.g. |n^e(\chi)| in eq. (35) of Uphoff et al. (2023))
  // is incorporated. This explains the factor 2 (doubledSurfaceArea)
  //
  // Uphoff, C., May, D. A., & Gabriel, A. A. (2023). A discontinuous Galerkin method for
  // sequences of earthquakes and aseismic slip on multiple faults using unstructured curvilinear
  // grids. Geophysical Journal International, 233(1), 586-626.
  double doubledSurfaceArea{};
};

template <typename Cfg>
struct DREnergyOutput {
  Real<Cfg> slip[seissol::tensor::slipInterpolated<Cfg>::size()]{};
  // Sums over the whole simulation of increments that do not change sign, of the slip rate
  // magnitude and of the frictional power; in single precision, the small ones after the slip has
  // stopped would fall below the rounding of the sums and be lost, which biases the seismic moment
  // and the potency low
  double accumulatedSlip[seissol::dr::misc::NumPaddedPoints<Cfg>]{};
  double frictionalEnergy[seissol::dr::misc::NumPaddedPoints<Cfg>]{};
  Real<Cfg> timeSinceSlipRateBelowThreshold[seissol::dr::misc::NumPaddedPoints<Cfg>]{};

  static std::vector<seissol::io::datatype::StructDatatype::MemberInfo> datatypeLayout() {
    return {
        seissol::io::datatype::StructDatatype::MemberInfo{
            "slip",
            offsetof(DREnergyOutput, slip),
            seissol::io::datatype::inferDatatype<decltype(slip)>()},
        seissol::io::datatype::StructDatatype::MemberInfo{
            "accumulatedSlip",
            offsetof(DREnergyOutput, accumulatedSlip),
            seissol::io::datatype::inferDatatype<decltype(accumulatedSlip)>()},
        seissol::io::datatype::StructDatatype::MemberInfo{
            "frictionalEnergy",
            offsetof(DREnergyOutput, frictionalEnergy),
            seissol::io::datatype::inferDatatype<decltype(frictionalEnergy)>()},
        seissol::io::datatype::StructDatatype::MemberInfo{
            "timeSinceSlipRateBelowThreshold",
            offsetof(DREnergyOutput, timeSinceSlipRateBelowThreshold),
            seissol::io::datatype::inferDatatype<decltype(timeSinceSlipRateBelowThreshold)>()},
    };
  }
};

} // namespace seissol

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_TYPEDEFS_H_
