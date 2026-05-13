// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_MISC_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_MISC_H_

#include "Common/Constants.h"
#include "Common/Marker.h"
#include "GeneratedCode/init.h"
#include "Geometry/MeshDefinition.h"
#include "Kernels/Precision.h"
#include "Solver/MultipleSimulations.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <tuple>
#include <type_traits>

namespace seissol::dr::misc {

/**
 * Stores the different types of friction laws
 * The values resemble the identifiers used in the old fortran implementation.
 * And they are still used in the parameter file like that.
 */
enum class FrictionLawType : uint32_t {
  NoFault = 0,
  LinearSlipWeakeningLegacy = 2,
  LinearSlipWeakening = 16,
  LinearSlipWeakeningBimaterial = 6,
  LinearSlipWeakeningTPApprox = 1058,
  RateAndStateAgingLaw = 3,
  RateAndStateSlipLaw = 4,
  RateAndStateFastVelocityWeakening = 103,
  ImposedSlipRatesYoffe = 33,
  ImposedSlipRatesGaussian = 34,
  ImposedSlipRatesDelta = 35,
  RateAndStateSevereVelocityWeakening = 7,
  RateAndStateAgingNucleation = 101,
};

// TODO: this can be moved to yateto headers
template <typename Tensor, int Dim>
constexpr uint32_t dimSize() noexcept {
  return Tensor::Stop[Dim] - Tensor::Start[Dim];
}

template <typename Tensor>
constexpr uint32_t leadDim() noexcept {
  return dimSize<Tensor, 0>();
}

/**
 * Number of gauss points padded to match the vector register length.
 */
static constexpr inline uint32_t NumPaddedPoints =
    multisim::MultisimEnabled
        ? dimSize<init::QInterpolated, 0>() * dimSize<init::QInterpolated, 1>()
        : leadDim<init::QInterpolated>();
static constexpr inline uint32_t NumPaddedPointsSingleSim =
    dimSize<init::QInterpolated, multisim::BasisFunctionDimension>();
static constexpr inline uint32_t NumQuantities =
    misc::dimSize<init::QInterpolated, multisim::BasisFunctionDimension + 1>();

/*
 * Time integration point count
 */

static constexpr inline uint32_t TimeSteps = ConvergenceOrder;

/**
 * Constants for Thermal Pressurization
 */
static constexpr uint32_t NumTpGridPoints = 60;
static constexpr double TpLogDz = 0.3;
static constexpr double TpMaxWaveNumber = 10.0;

/**
 * Number of gauss points on an element surface.
 */
static constexpr uint32_t NumBoundaryGaussPoints =
    init::QInterpolated::Shape[multisim::BasisFunctionDimension];

template <std::size_t I, typename F, typename TupleT>
constexpr F forEachElement(F&& functor, TupleT&& tuple) {
  // TODO: maybe forward here somehow?
  functor(std::get<I>(tuple), I);
  if constexpr (I + 1 < std::tuple_size_v<std::remove_reference_t<TupleT>>) {
    return forEachElement<I + 1>(std::forward<F>(functor), std::forward<TupleT>(tuple));
  }
  return std::forward<F>(functor);
}

template <typename TupleT, typename F>
constexpr F forEach(TupleT&& tuple, F&& functor) {
  return forEachElement<0>(std::forward<F>(functor), std::forward<TupleT>(tuple));
}

/**
 * Compute base^exp
 * @param base
 * @return
 */
#pragma omp declare simd
template <int32_t Exp, typename T>
SEISSOL_HOSTDEVICE constexpr auto power(T base) -> T {
  if constexpr (Exp == 0) {
    return 1;
  }
  if constexpr (Exp == 1) {
    return base;
  }
  if constexpr (Exp == 2) {
    return base * base;
  }

  constexpr std::int32_t ILogExp = [&]() constexpr {
    auto val = Exp;
    for (std::int32_t i = 0; i < 64; ++i) {
      val /= 2;
      if (val == 0) {
        // (Exp == 0 is handled above)
        return i;
      }
    }
    return 64;
  }();

  T res = base;

#ifdef ACL_DEVICE
#pragma unroll
#endif
  for (std::int32_t i = ILogExp - 1; i >= 0; --i) {
    res *= res;
    if ((Exp & (1 << i)) != 0) {
      res *= base;
    }
  }

  return res;
}

#pragma omp declare simd
template <typename T>
SEISSOL_HOSTDEVICE constexpr T square(T t) {
  return t * t;
}

/**
 * Computes a squared sum of an N-dimensional vector
 * @return magnitude of the vector
 */
#pragma omp declare simd
template <typename T, typename... Tn>
SEISSOL_HOSTDEVICE constexpr T square(T t1, Tn... tn) {
  return square(t1) + square(tn...);
}

/**
 * Computes the magnitude of an N-dimensional vector
 * @return magnitude of the vector
 */
#pragma omp declare simd
template <typename T, typename... Tn>
SEISSOL_HOSTDEVICE constexpr T magnitude(T t1, Tn... tn) {
  static_assert((std::is_same_v<T, Tn> && ...), "All types need to be equal.");
  if constexpr (sizeof...(Tn) == 1) {
    return std::hypot(t1, tn...);
  }
  return std::sqrt(square(t1, tn...));
}

#pragma omp declare simd
template <typename T>
SEISSOL_HOSTDEVICE constexpr T clamp(T value, T minval, T maxval) {
#ifdef __HIP__
  return std::max(minval, std::min(maxval, value));
#else
  return std::clamp(value, minval, maxval);
#endif
}

/**
 * Create strike and dip unit vectors give a fault normal vector
 * Note: equations are explained in documentation -> left-lateral-right-lateral-normal-reverse
 * @param normal
 * @param strike
 * @param dip
 */
void computeStrikeAndDipVectors(const VrtxCoords normal, VrtxCoords strike, VrtxCoords dip);

std::string frictionLawName(seissol::dr::misc::FrictionLawType type);

// NOLINTBEGIN (-cppcoreguidelines-use-enum-class)

namespace quantity_indices {
/**
 * Defines the indices under which one can find a specific quantity.
 * U, V, W: Velocities in x, y, z direction.
 * N, T1, T2: traction in normal and fault aligned directions.
 * XX, YY, ZZ, XY, YZ, XZ: Stress in cartesian coordinates
 * Use as:
 * ```
 * using namepace dr::misc::quantity_indices;
 * real quantities[9];
 * real normalStress = quantities[N];
 * ```
 * */
enum QuantityIndices : uint32_t {
  U = 6,
  V = 7,
  W = 8,
  N = 0,
  T1 = 3,
  T2 = 5,
  XX = 0,
  YY = 1,
  ZZ = 2,
  XY = 3,
  YZ = 4,
  XZ = 5,
};
} // namespace quantity_indices
} // namespace seissol::dr::misc

// NOLINTEND ()

namespace seissol::initializer::parameters {
struct DRParameters;
} // namespace seissol::initializer::parameters

namespace seissol::dr {
// compile-time parameter; rather arbitrary (and just large enough for most cases). It's there to
// avoid us allocating dynamic arrays in the parameters.
constexpr std::size_t MaxNucleations = 16;

/**
 * Friction law parameters, as used in the kernels.
 * For separation of concerns (and using the `real` datatype), prefer this one
 * to the one in the Initializer/Parameters.
 */
struct FrictionLawParameters {
  real healingThreshold{-1.0};
  real tpProxyExponent{0.0};
  real rsF0{0.0};
  real rsB{0.0};
  real rsSr0{0.0};
  real rsInitialSlipRate1{0.0};
  real rsInitialSlipRate2{0.0};
  real muW{0.0};
  real thermalDiffusivity{0.0};
  real heatCapacity{0.0};
  real undrainedTPResponse{0.0};
  real initialTemperature{0.0};
  real initialPressure{0.0};
  // Prakash-Clifton regularization parameter
  real vStar{0.0};
  real prakashLength{0.0};
  real terminatorSlipRateThreshold{0.0};
  real etaDamp{1.0};
  real etaDampEnd{std::numeric_limits<real>::infinity()};
  std::array<real, MaxNucleations> t0{};
  std::array<real, MaxNucleations> s0{};
  std::uint32_t nucleationCount{0};
  std::uint32_t rsMaxNumberSlipRateUpdates{60};
  std::uint32_t rsNumberStateVariableUpdates{10};
  real rsSlipRateTolerance{1e-8};
  real rsStateTolerance{1e-8};
  bool isFrictionEnergyRequired{false};
  bool isCheckAbortCriteraEnabled{false};
  bool energiesFromAcrossFaultVelocities{false};

  FrictionLawParameters() = default;
  explicit FrictionLawParameters(const seissol::initializer::parameters::DRParameters& parameters);
};
} // namespace seissol::dr

#endif // SEISSOL_SRC_DYNAMICRUPTURE_MISC_H_
