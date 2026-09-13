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
#include "Common/Marker.h"
#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "Kernels/Precision.h"

#include <cmath>
#include <limits>

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

template <Executor Executor>
struct ImposedState;

/// Whether the impedance of a face varies from one of its nodes to the next.
///
/// It does where the moduli a wave sees follow the state rather than the
/// material: the impedance is then a property of a point and an instant, and
/// there is nothing about it to keep between timesteps.
inline constexpr bool NodalImpedance = model::MaterialT::Type == model::MaterialType::Damage;

/// What the impedance of a node needs from the material of the two sides.
///
/// Everything the effective moduli need except the state of the node itself.
/// The density and the undamaged moduli could be read back out of the wave
/// speeds of the face, but only if those mean the undamaged ones -- they are
/// kept here instead, where what they mean is written down. Not stored at all
/// where the impedance belongs to the face.
struct NodalImpedanceParameters {
  double rhoPlus{};
  double lambda0Plus{};
  double mu0Plus{};
  double gammaRPlus{};
  double xi0Plus{};
  double rhoMinus{};
  double lambda0Minus{};
  double mu0Minus{};
  double gammaRMinus{};
  double xi0Minus{};
};

/// Fills what the impedance of a node needs from the material of the two
/// sides of a face. Nothing to fill where the impedance belongs to the face.
///
/// A template so that the members it reads are looked up when a material is
/// put in, and not before: no material whose impedance belongs to its face
/// has them, and there is no build in which both kinds are compiled.
template <typename MaterialT>
void setNodalImpedanceParameters([[maybe_unused]] const MaterialT& plus,
                                 [[maybe_unused]] const MaterialT& minus,
                                 [[maybe_unused]] NodalImpedanceParameters& parameters) {
  if constexpr (NodalImpedance) {
    parameters.rhoPlus = plus.rho;
    parameters.lambda0Plus = plus.lambda0;
    parameters.mu0Plus = plus.mu0;
    parameters.gammaRPlus = plus.gammaR;
    parameters.xi0Plus = plus.xi0;
    parameters.rhoMinus = minus.rho;
    parameters.lambda0Minus = minus.lambda0;
    parameters.mu0Minus = minus.mu0;
    parameters.gammaRMinus = minus.gammaR;
    parameters.xi0Minus = minus.xi0;
  }
}

/**
 * The impedances of a face's nodes, where they belong to the nodes.
 *
 * Empty for every material whose moduli come from the material, and then free:
 * the fault stresses derive from it, and an empty base costs no space. Where it
 * is not empty it lives exactly as the stresses do -- an array over the nodes
 * on the host, one value per thread on the device, made in the precomputation
 * and gone with it. Deriving rather than standing beside them is what keeps
 * every signature that passes the stresses unchanged, which is every hook of
 * every friction law.
 */
template <Executor Executor, bool Nodal>
struct FaultImpedancesImpl {};

template <>
struct FaultImpedancesImpl<Executor::Host, true> {
  alignas(Alignment) real etaS[misc::NumPaddedPoints]{};
  alignas(Alignment) real invEtaS[misc::NumPaddedPoints]{};
  alignas(Alignment) real invZs[misc::NumPaddedPoints]{};
  alignas(Alignment) real invZp[misc::NumPaddedPoints]{};
  alignas(Alignment) real invZsNeig[misc::NumPaddedPoints]{};
  alignas(Alignment) real invZpNeig[misc::NumPaddedPoints]{};
};

template <>
struct FaultImpedancesImpl<Executor::Device, true> {
  real etaS{};
  real invEtaS{};
  real invZs{};
  real invZp{};
  real invZsNeig{};
  real invZpNeig{};
};

template <Executor Executor>
using FaultImpedances = FaultImpedancesImpl<Executor, NodalImpedance>;

/// The impedance at one node, which is what the device keeps per thread. The
/// same type, so that a node's impedance and a face's cannot drift into two
/// definitions of the same six numbers.
using NodalImpedanceT = FaultImpedancesImpl<Executor::Device, true>;

/// The impedance a wave sees at one node of a face.
///
/// The moduli follow the state here, so this is a property of a point and an
/// instant rather than of the material. It is the secant linearisation the
/// volume carries, not the exact tangent: the tangent has terms in the outer
/// product of the strain with itself, and its acoustic tensor would depend on
/// the direction of the normal relative to the principal strain axes. The
/// secant is what the volume's wave speed uses, and a face that carried a
/// different material than the cells beside it is worse than a face that
/// carries an approximate one.
///
/// The solid branch alone, without the breakage blend, for the same reason.
SEISSOL_HOSTDEVICE inline NodalImpedanceT nodalImpedance(const NodalImpedanceParameters& params,
                                                         real alphaPlus,
                                                         real xiPlus,
                                                         real alphaMinus,
                                                         real xiMinus) {
  const auto shear = [](double mu0, double gammaR, double xi0, real alpha, real xi) {
    // 2 mu_eff, as the volume forms it
    return static_cast<real>(2.0 * mu0 - 2.0 * gammaR * xi0 * alpha - gammaR * alpha * xi);
  };

  const auto twoMuPlus =
      shear(params.mu0Plus, params.gammaRPlus, params.xi0Plus, alphaPlus, xiPlus);
  const auto twoMuMinus =
      shear(params.mu0Minus, params.gammaRMinus, params.xi0Minus, alphaMinus, xiMinus);

  const auto zp = std::sqrt(static_cast<real>(params.rhoPlus) *
                            (static_cast<real>(params.lambda0Plus) + twoMuPlus));
  const auto zpNeig = std::sqrt(static_cast<real>(params.rhoMinus) *
                                (static_cast<real>(params.lambda0Minus) + twoMuMinus));
  const auto zs = std::sqrt(static_cast<real>(params.rhoPlus) * static_cast<real>(0.5) * twoMuPlus);
  const auto zsNeig =
      std::sqrt(static_cast<real>(params.rhoMinus) * static_cast<real>(0.5) * twoMuMinus);

  NodalImpedanceT impedance{};
  impedance.invZp = static_cast<real>(1.0) / zp;
  impedance.invZpNeig = static_cast<real>(1.0) / zpNeig;
  impedance.invZs = static_cast<real>(1.0) / zs;
  impedance.invZsNeig = static_cast<real>(1.0) / zsNeig;
  impedance.invEtaS = impedance.invZs + impedance.invZsNeig;
  impedance.etaS = static_cast<real>(1.0) / impedance.invEtaS;
  return impedance;
}

/// The strain invariant ratio at a node, from the six Voigt components of the
/// strain there. An invariant, so the rotation into the face frame does not
/// enter and a face may form it from the rotated strain directly.
SEISSOL_HOSTDEVICE inline real
    strainRatio(real exx, real eyy, real ezz, real exy, real eyz, real exz) {
  const auto i1 = exx + eyy + ezz;
  const auto i2 = exx * exx + eyy * eyy + ezz * ezz +
                  static_cast<real>(2.0) * (exy * exy + eyz * eyz + exz * exz);
  const auto floor = std::numeric_limits<real>::epsilon() * std::numeric_limits<real>::epsilon();
  return i2 > floor ? i1 / std::sqrt(i2) : static_cast<real>(0.0);
}

/**
 * Struct that contains all input stresses
 * normalStress in direction of the face normal, traction1, traction2 in the direction of the
 * respective tangential vectors
 */
template <>
struct FaultStresses<Executor::Host> : FaultImpedances<Executor::Host> {
  alignas(Alignment) real normalStress[misc::NumPaddedPoints]{};
  alignas(Alignment) real traction1[misc::NumPaddedPoints]{};
  alignas(Alignment) real traction2[misc::NumPaddedPoints]{};
  alignas(Alignment) real fluidPressure[misc::NumPaddedPoints]{};
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
  alignas(Alignment) real normalStress[misc::NumPaddedPoints]{};
  alignas(Alignment) real traction1[misc::NumPaddedPoints]{};
  alignas(Alignment) real traction2[misc::NumPaddedPoints]{};
};

/**
 * Accumulator for the imposed state. Used after every internal timestep.
 */
template <>
struct ImposedState<Executor::Host> {
  alignas(Alignment) real plus[misc::NumQuantities][misc::NumPaddedPoints]{};
  alignas(Alignment) real minus[misc::NumQuantities][misc::NumPaddedPoints]{};
};

/**
 * Struct that contains all input stresses
 * normalStress in direction of the face normal, traction1, traction2 in the direction of the
 * respective tangential vectors
 */
template <>
struct FaultStresses<Executor::Device> : FaultImpedances<Executor::Device> {
  real normalStress{};
  real traction1{};
  real traction2{};
  real fluidPressure{};
};

/**
 * Struct that contains all traction results
 * normalStress in direction of the face normal, traction1, traction2 in the direction of the
 * respective tangential vectors. See TractionResults<Executor::Host> for the semantics of
 * normalStress.
 */
template <>
struct TractionResults<Executor::Device> {
  real normalStress{};
  real traction1{};
  real traction2{};
};

/**
 * Accumulator for the imposed state. Used after every internal timestep.
 */
template <>
struct ImposedState<Executor::Device> {
  real plus[misc::NumQuantities]{};
  real minus[misc::NumQuantities]{};
};

} // namespace seissol::dr

#endif // SEISSOL_SRC_DYNAMICRUPTURE_TYPEDEFS_H_
