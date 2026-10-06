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
#include "Kernels/Precision.h"
#include "Model/OperatorLayout.h"

#include <cstddef>

namespace seissol::dr {

/**
 * One scalar of the Riemann problem, read at a point of the fault.
 *
 * A material that varies along a face gives a different value at every
 * quadrature point of it; one that does not gives the same value everywhere and
 * has no reason to store it more than once. Both answer the same question --
 * what is this scalar at point `index` -- so a friction law is written once and
 * the storage follows the material. Where the value is one number the index
 * never reaches the array and the compiler hoists the read out of the point
 * loop again, which is what it did when the scalar was a plain member.
 */
template <bool Pointwise, typename T = real>
class PointScalar {
  public:
  static constexpr std::size_t Count = Pointwise ? misc::NumPaddedPoints : 1;

#pragma omp declare simd
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr T operator()(std::size_t index) const {
    return values[Pointwise ? index : 0];
  }

  /// The same value at every point of the face.
  SEISSOL_HOSTDEVICE constexpr void fill(T value) {
    for (std::size_t point = 0; point < Count; ++point) {
      values[point] = value;
    }
  }

  SEISSOL_HOSTDEVICE constexpr void set(std::size_t index, T value) {
    values[Pointwise ? index : 0] = value;
  }

  private:
  T values[Count]{};
};

/**
 * Stores the P and S wave impedances for an element and its neighbor as well as the eta values from
 * Carsten Uphoff's dissertation equation (4.51)
 */
template <bool Pointwise>
struct ImpedancesAndEtaOf {
  PointScalar<Pointwise> zp;
  PointScalar<Pointwise> zs;
  PointScalar<Pointwise> zpNeig;
  PointScalar<Pointwise> zsNeig;
  PointScalar<Pointwise> etaP;
  PointScalar<Pointwise> etaS;
  PointScalar<Pointwise> invEtaS;
  PointScalar<Pointwise> invZp;
  PointScalar<Pointwise> invZs;
  PointScalar<Pointwise> invZpNeig;
  PointScalar<Pointwise> invZsNeig;
};

/// Whether the fault reads its scalars per point, which it does exactly when the
/// material is allowed to vary within a cell.
constexpr bool PointwiseImpedances = NodalMaterial;

using ImpedancesAndEta = ImpedancesAndEtaOf<PointwiseImpedances>;

/// How many points of a face carry their own Riemann problem.
constexpr std::size_t ImpedancePoints = PointScalar<PointwiseImpedances>::Count;

/**
 * The density and the wave speeds on one side of a fault, read at a point of it. What the fault
 * derives from the material beyond its Riemann problem -- the modulus that turns slip into
 * seismic moment, the stress components outside the Riemann problem the receivers reconstruct --
 * reads them here, so that it sees the same material at a point as the impedances do. Stored in
 * double in every build, as the members of the material are, where the impedances are real. Where
 * the material varies inside a cell, the material at a point comes out of the generated
 * projection of the samples, which computes in real, so in a single precision build these values
 * carry no more than single precision there.
 */
template <bool Pointwise>
struct WaveSpeedsOf {
  PointScalar<Pointwise, double> density;
  PointScalar<Pointwise, double> pWaveVelocity;
  PointScalar<Pointwise, double> sWaveVelocity;
};

using WaveSpeeds = WaveSpeedsOf<PointwiseImpedances>;

namespace internal {
/// A count of reals rounded up to the vector width.
constexpr std::size_t paddedToVector(std::size_t count) {
  constexpr std::size_t Width = Alignment / sizeof(real);
  return (count + Width - 1) / Width * Width;
}
} // namespace internal

/**
 * How a fault face stores the lift of one side: the operator that takes the imposed state of that
 * side, given in the coordinates of the face at its quadrature points, into the side's cell.
 *
 * Where the material does not vary inside a cell, that is one matrix, fluxSolver, with the scale
 * of the side and the rotation back to global coordinates folded in. Where it does, the operator
 * differs from point to point; the face then stores its rotation once and, for every point, the
 * scalars the coefficient matrix of the fault normal is linear in, scaled. Each block starts at a
 * multiple of the vector width.
 */
struct FaultFluxLayout {
  /// Where the rotation starts, in reals.
  static constexpr std::size_t RotationOffset = 0;

  /// How many reals the rotation takes: one for the face, or one per point where it may be curved.
  static constexpr std::size_t RotationSize =
      Curvilinear ? tensor::TPoints::size() : tensor::T::size();

  /// Where the scalars of the first coefficient start.
  static constexpr std::size_t CoefficientsOffset = internal::paddedToVector(RotationSize);

  /// How far apart the scalars of two coefficients are.
  static constexpr std::size_t CoefficientStride =
      internal::paddedToVector(misc::NumBoundaryGaussPoints);

  /// Where the scalars of one coefficient start, one per quadrature point of the face.
  static constexpr std::size_t coefficientOffset(std::size_t coefficient) {
    return CoefficientsOffset + coefficient * CoefficientStride;
  }

  /// How many reals the lift of one side takes.
  static constexpr std::size_t Size =
      NodalFaultFlux ? CoefficientsOffset + FaultFluxCoefficientCount * CoefficientStride
                     : tensor::fluxSolver::size();
};

/**
 * One matrix of the Riemann problem, read at a point of the fault. The
 * counterpart of PointScalar for the quantities that are not scalars; a point
 * gets a contiguous matrix so that a caller can keep reading it as one.
 */
template <bool Pointwise, std::size_t Size>
class PointMatrix {
  public:
  static constexpr std::size_t Count = Pointwise ? misc::NumPaddedPoints : 1;

  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr const real* at(std::size_t index) const {
    return values[Pointwise ? index : 0];
  }

  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr real* at(std::size_t index) {
    return values[Pointwise ? index : 0];
  }

  private:
  alignas(Alignment) real values[Count][Size]{};
};

/**
 * Stores the matrices of the Riemann problem of a fault face, for an element and its neighbor.
 * The admittances, eta and the lateral stress map generalize equation (4.51) from Carsten's
 * thesis to a Riemann problem that couples the traction components; they are filled only for the
 * materials that have one, anisotropic and poroelastic ones. The traction averaging matrices are
 * filled for every material.
 */
template <bool Pointwise>
struct ImpedanceMatricesOf {
  PointMatrix<Pointwise, tensor::Zplus::size()> impedance;
  PointMatrix<Pointwise, tensor::Zminus::size()> impedanceNeig;
  PointMatrix<Pointwise, tensor::eta::size()> eta;
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
  PointMatrix<Pointwise, 3 * tensor::Zminus::Shape[0]> lateralStress;
  /**
   * The traction averaging matrices b+ and b-: the weights with which the state of either side
   * enters the traction of the interface, tau = b+^T q+ + b-^T q-, stored in the layout of the
   * tractionPlusMatrix and tractionMinusMatrix tensors. The static frictional work of the energy
   * output contracts them for every material; the friction energy contracts them only for an
   * anisotropic one, where the normal traction enters the shear traction, and takes the same
   * weights from the scalar impedances otherwise. They follow from the impedances of a point
   * exactly as eta does, so they vary along the face wherever the impedances do.
   */
  PointMatrix<Pointwise, tensor::tractionPlusMatrix::size()> tractionPlus;
  PointMatrix<Pointwise, tensor::tractionMinusMatrix::size()> tractionMinus;
};

using ImpedanceMatrices = ImpedanceMatricesOf<PointwiseImpedances>;

template <Executor Executor>
struct FaultStresses;

template <Executor Executor>
struct TractionResults;

template <Executor Executor>
struct ImposedState;

/**
 * Struct that contains all input stresses
 * normalStress in direction of the face normal, traction1, traction2 in the direction of the
 * respective tangential vectors
 */
template <>
struct FaultStresses<Executor::Host> {
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
struct FaultStresses<Executor::Device> {
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
