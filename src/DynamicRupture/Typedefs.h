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

#include <array>
#include <cmath>
#include <cstddef>
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
  real rhoPlus{};
  real lambda0Plus{};
  real mu0Plus{};
  real gammaRPlus{};
  real xi0Plus{};
  real rhoMinus{};
  real lambda0Minus{};
  real mu0Minus{};
  real gammaRMinus{};
  real xi0Minus{};
  /// The strain each side carries before the first timestep, rotated into the
  /// frame of this face. The tangent of the stress turns on the direction of
  /// the total strain, and at a fault the background dominates it.
  real epsInitPlus[6]{};
  real epsInitMinus[6]{};
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
                                 [[maybe_unused]] const std::array<double, 36>& rotation,
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

    // A strain in tensor components rotates with the matrix a stress rotates
    // with; only the engineering convention, which doubles the shear, would
    // need the other one.
    const double plusStrain[6] = {plus.epsInitXX,
                                  plus.epsInitYY,
                                  plus.epsInitZZ,
                                  plus.epsInitXY,
                                  plus.epsInitYZ,
                                  plus.epsInitXZ};
    const double minusStrain[6] = {minus.epsInitXX,
                                   minus.epsInitYY,
                                   minus.epsInitZZ,
                                   minus.epsInitXY,
                                   minus.epsInitYZ,
                                   minus.epsInitXZ};
    for (std::size_t row = 0; row < 6; ++row) {
      double rotatedPlus = 0.0;
      double rotatedMinus = 0.0;
      for (std::size_t col = 0; col < 6; ++col) {
        rotatedPlus += rotation[row * 6 + col] * plusStrain[col];
        rotatedMinus += rotation[row * 6 + col] * minusStrain[col];
      }
      parameters.epsInitPlus[row] = static_cast<real>(rotatedPlus);
      parameters.epsInitMinus[row] = static_cast<real>(rotatedMinus);
    }
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
  alignas(Alignment) real admittance[9][misc::NumPaddedPoints]{};
  alignas(Alignment) real admittanceNeig[9][misc::NumPaddedPoints]{};
  alignas(Alignment) real eta[9][misc::NumPaddedPoints]{};
  alignas(Alignment) real lateralStress[9][misc::NumPaddedPoints]{};
};

template <>
struct FaultImpedancesImpl<Executor::Device, true> {
  real admittance[9]{};
  real admittanceNeig[9]{};
  real eta[9]{};
  real lateralStress[9]{};
};

template <Executor Executor>
using FaultImpedances = FaultImpedancesImpl<Executor, NodalImpedance>;

/// The impedance at one node, which is what the device keeps per thread. The
/// same type, so that a node's impedance and a face's cannot drift into two
/// definitions of the same three matrices.
using NodalImpedanceT = FaultImpedancesImpl<Executor::Device, true>;

/// Eigenvalues of a symmetric 3x3, from the trigonometric solution of its
/// characteristic cubic. No iteration, and no branch beyond the one that says
/// the matrix is already diagonal.
SEISSOL_HOSTDEVICE inline void eigenvaluesSymmetric(const real matrix[9], real values[3]) {
  const auto offDiagonal = matrix[1] * matrix[1] + matrix[2] * matrix[2] + matrix[5] * matrix[5];
  if (offDiagonal <= static_cast<real>(0.0)) {
    values[0] = matrix[0];
    values[1] = matrix[4];
    values[2] = matrix[8];
    return;
  }
  const auto mean = (matrix[0] + matrix[4] + matrix[8]) / static_cast<real>(3.0);
  const auto spread =
      (matrix[0] - mean) * (matrix[0] - mean) + (matrix[4] - mean) * (matrix[4] - mean) +
      (matrix[8] - mean) * (matrix[8] - mean) + static_cast<real>(2.0) * offDiagonal;
  const auto scale = std::sqrt(spread / static_cast<real>(6.0));
  real shifted[9];
  for (std::size_t k = 0; k < 9; ++k) {
    shifted[k] = matrix[k] / scale;
  }
  shifted[0] -= mean / scale;
  shifted[4] -= mean / scale;
  shifted[8] -= mean / scale;
  const auto determinant = shifted[0] * (shifted[4] * shifted[8] - shifted[5] * shifted[7]) -
                           shifted[1] * (shifted[3] * shifted[8] - shifted[5] * shifted[6]) +
                           shifted[2] * (shifted[3] * shifted[7] - shifted[4] * shifted[6]);
  auto argument = determinant / static_cast<real>(2.0);
  argument = argument < static_cast<real>(-1.0)
                 ? static_cast<real>(-1.0)
                 : (argument > static_cast<real>(1.0) ? static_cast<real>(1.0) : argument);
  const auto angle = std::acos(argument) / static_cast<real>(3.0);
  constexpr auto TwoThirdsPi = static_cast<real>(2.0943951023931953);
  values[0] = mean + static_cast<real>(2.0) * scale * std::cos(angle);
  values[2] = mean + static_cast<real>(2.0) * scale * std::cos(angle + TwoThirdsPi);
  values[1] = static_cast<real>(3.0) * mean - values[0] - values[2];
}

/// The inverse square root of a symmetric positive definite 3x3.
///
/// Through Cayley-Hamilton on both halves: the square root is a combination of
/// the matrix, its square and the identity once the invariants of the root are
/// known, and inverting it is that statement again. No eigenvector is formed
/// and no inverse is taken, which is what makes it exact where two of the
/// three eigenvalues coincide -- at a fault, the two shear waves.
SEISSOL_HOSTDEVICE inline void inverseSqrtSymmetric(const real matrix[9], real result[9]) {
  const auto product = [](const real* a, const real* b, real* out) {
    for (std::size_t row = 0; row < 3; ++row) {
      for (std::size_t col = 0; col < 3; ++col) {
        real sum = static_cast<real>(0.0);
        for (std::size_t k = 0; k < 3; ++k) {
          sum += a[3 * row + k] * b[3 * k + col];
        }
        out[3 * row + col] = sum;
      }
    }
  };
  real values[3];
  eigenvaluesSymmetric(matrix, values);
  const auto first = std::sqrt(values[0]) + std::sqrt(values[1]) + std::sqrt(values[2]);
  const auto third = std::sqrt(values[0] * values[1] * values[2]);
  const auto trace = matrix[0] + matrix[4] + matrix[8];
  const auto second = static_cast<real>(0.5) * (first * first - trace);

  real square[9];
  product(matrix, matrix, square);
  const auto scale = first * second - third;
  real root[9];
  for (std::size_t k = 0; k < 9; ++k) {
    root[k] = (-square[k] + (first * first - second) * matrix[k]) / scale;
  }
  root[0] += first * third / scale;
  root[4] += first * third / scale;
  root[8] += first * third / scale;

  real rootSquare[9];
  product(root, root, rootSquare);
  for (std::size_t k = 0; k < 9; ++k) {
    result[k] = (rootSquare[k] - first * root[k]) / third;
  }
  result[0] += second / third;
  result[4] += second / third;
  result[8] += second / third;
}

/// The tangent of the stress at one node, in Voigt, as the 6 by 6 that a
/// strain in tensor components is contracted with.
///
/// With n the strain normalised in the Frobenius norm, xi = tr(n) and
/// g = gammaR alpha,
///
///   C = lambda0 d(x)d + (2 mu0 - 2 g xi0 - g xi) Isym
///         - g (d(x)n + n(x)d) + g xi n(x)n,
///
/// whose last two groups no pair of Lame parameters expresses. The solid
/// branch alone, without the breakage blend.
SEISSOL_HOSTDEVICE inline void damageTangent(real lambda0,
                                             real mu0,
                                             real gammaR,
                                             real xi0,
                                             real alpha,
                                             const real strain[6],
                                             real tangent[36]) {
  constexpr real Weight[6] = {static_cast<real>(1.0),
                              static_cast<real>(1.0),
                              static_cast<real>(1.0),
                              static_cast<real>(2.0),
                              static_cast<real>(2.0),
                              static_cast<real>(2.0)};
  const auto i1 = strain[0] + strain[1] + strain[2];
  real i2 = static_cast<real>(0.0);
  for (std::size_t k = 0; k < 6; ++k) {
    i2 += Weight[k] * strain[k] * strain[k];
  }
  const auto floor = std::numeric_limits<real>::epsilon() * std::numeric_limits<real>::epsilon();
  const auto deformed = i2 > floor;
  const auto root = std::sqrt(deformed ? i2 : floor);
  // A node at rest has no direction to give, and there the material is the
  // undamaged one whatever the damage says.
  const auto xi = deformed ? i1 / root : static_cast<real>(0.0);
  const auto coupling = deformed ? gammaR * alpha : static_cast<real>(0.0);
  const auto shear =
      static_cast<real>(2.0) * mu0 - static_cast<real>(2.0) * coupling * xi0 - coupling * xi;

  constexpr real Delta[6] = {static_cast<real>(1.0),
                             static_cast<real>(1.0),
                             static_cast<real>(1.0),
                             static_cast<real>(0.0),
                             static_cast<real>(0.0),
                             static_cast<real>(0.0)};
  for (std::size_t a = 0; a < 6; ++a) {
    const auto na = deformed ? strain[a] / root : static_cast<real>(0.0);
    for (std::size_t b = 0; b < 6; ++b) {
      const auto nb = deformed ? strain[b] / root : static_cast<real>(0.0);
      auto value = lambda0 * Delta[a] * Delta[b] - coupling * (Delta[a] * nb + na * Delta[b]) +
                   coupling * xi * na * nb;
      if (a == b) {
        value += shear * (a < 3 ? static_cast<real>(1.0) : static_cast<real>(0.5));
      }
      tangent[6 * a + b] = value;
    }
  }
}

/// The admittance a wave sees at one node of a face.
///
/// The acoustic tensor of the tangent in the face normal, which the face frame
/// puts on the first axis, is Gamma_ik = C_i1k1; the admittance is
/// (rho Gamma)^-1/2, which is where a scalar would carry 1 / (rho c).
SEISSOL_HOSTDEVICE inline void
    admittanceFromTangent(real rho, const real tangent[36], real admittance[9]) {
  // Voigt index of the pair (i, normal): nn, nt1, nt2. Read as plain tensor
  // components -- a Voigt weight belongs to a contraction with a strain, and
  // this is not one.
  constexpr std::size_t Traction[3] = {0, 3, 5};
  real gamma[9];
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t k = 0; k < 3; ++k) {
      gamma[3 * i + k] = rho * tangent[6 * Traction[i] + Traction[k]];
    }
  }
  inverseSqrtSymmetric(gamma, admittance);
}

/// The inverse of a 3x3, by cofactors.
SEISSOL_HOSTDEVICE inline void inverse3(const real matrix[9], real result[9]) {
  const auto determinant = matrix[0] * (matrix[4] * matrix[8] - matrix[5] * matrix[7]) -
                           matrix[1] * (matrix[3] * matrix[8] - matrix[5] * matrix[6]) +
                           matrix[2] * (matrix[3] * matrix[7] - matrix[4] * matrix[6]);
  result[0] = (matrix[4] * matrix[8] - matrix[5] * matrix[7]) / determinant;
  result[1] = (matrix[2] * matrix[7] - matrix[1] * matrix[8]) / determinant;
  result[2] = (matrix[1] * matrix[5] - matrix[2] * matrix[4]) / determinant;
  result[3] = (matrix[5] * matrix[6] - matrix[3] * matrix[8]) / determinant;
  result[4] = (matrix[0] * matrix[8] - matrix[2] * matrix[6]) / determinant;
  result[5] = (matrix[2] * matrix[3] - matrix[0] * matrix[5]) / determinant;
  result[6] = (matrix[3] * matrix[7] - matrix[4] * matrix[6]) / determinant;
  result[7] = (matrix[1] * matrix[6] - matrix[0] * matrix[7]) / determinant;
  result[8] = (matrix[0] * matrix[4] - matrix[1] * matrix[3]) / determinant;
}

/// The stress a face carries outside its Riemann problem, as the map from the
/// traction it solves for.
///
/// Across the waves of the face normal only the strains with a traction index
/// jump; the three tangential ones carry no flux through the face and stand
/// still. So with T the traction indices and L the lateral ones,
///
///   d sigma_L = C_LT w d eps_T   and   d sigma_T = C_TT w d eps_T,
///
/// and the map is (C_LT w) (C_TT w)^-1. For an undamaged solid it comes out
/// as lambda / (lambda + 2 mu) on the two normal rows and nothing else;
/// damage fills in the shear column, which no isotropic modulus reaches.
SEISSOL_HOSTDEVICE inline void lateralFromTangent(const real tangent[36], real lateralStress[9]) {
  constexpr std::size_t Traction[3] = {0, 3, 5};
  constexpr std::size_t Lateral[3] = {1, 2, 4};
  constexpr real Weight[3] = {
      static_cast<real>(1.0), static_cast<real>(2.0), static_cast<real>(2.0)};
  real tractionBlock[9];
  real lateralBlock[9];
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t k = 0; k < 3; ++k) {
      tractionBlock[3 * i + k] = tangent[6 * Traction[i] + Traction[k]] * Weight[k];
      lateralBlock[3 * i + k] = tangent[6 * Lateral[i] + Traction[k]] * Weight[k];
    }
  }
  real inverse[9];
  inverse3(tractionBlock, inverse);
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t k = 0; k < 3; ++k) {
      real sum = static_cast<real>(0.0);
      for (std::size_t j = 0; j < 3; ++j) {
        sum += lateralBlock[3 * i + j] * inverse[3 * j + k];
      }
      lateralStress[3 * i + k] = sum;
    }
  }
}

/// Everything one node of a face needs: the admittance of both sides, the eta
/// they make, and the map onto the stress outside the Riemann problem. All of
/// it out of the tangent, which is formed once per side and read three times.
SEISSOL_HOSTDEVICE inline void nodalAdmittance(const NodalImpedanceParameters& params,
                                               real alphaPlus,
                                               const real strainPlus[6],
                                               real alphaMinus,
                                               const real strainMinus[6],
                                               real admittancePlus[9],
                                               real admittanceMinus[9],
                                               real eta[9],
                                               real lateralStress[9]) {
  real totalPlus[6];
  real totalMinus[6];
  for (std::size_t k = 0; k < 6; ++k) {
    totalPlus[k] = strainPlus[k] + params.epsInitPlus[k];
    totalMinus[k] = strainMinus[k] + params.epsInitMinus[k];
  }

  real tangentPlus[36];
  real tangentMinus[36];
  damageTangent(params.lambda0Plus,
                params.mu0Plus,
                params.gammaRPlus,
                params.xi0Plus,
                alphaPlus,
                totalPlus,
                tangentPlus);
  damageTangent(params.lambda0Minus,
                params.mu0Minus,
                params.gammaRMinus,
                params.xi0Minus,
                alphaMinus,
                totalMinus,
                tangentMinus);
  admittanceFromTangent(params.rhoPlus, tangentPlus, admittancePlus);
  admittanceFromTangent(params.rhoMinus, tangentMinus, admittanceMinus);
  // The plus side, because that is where the fault output evaluates them.
  lateralFromTangent(tangentPlus, lateralStress);

  real sum[9];
  for (std::size_t k = 0; k < 9; ++k) {
    sum[k] = admittancePlus[k] + admittanceMinus[k];
  }
  inverse3(sum, eta);
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
