// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_FRICTIONLAWCONTEXT_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_FRICTIONLAWCONTEXT_H_

#include "Alignment.h"
#include "Core/Common/Executor.h"
#include "Core/Common/Marker.h"
#include "Core/Common/MultipleSimulations.h"
#include "Core/Common/Real.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolverCommon.h"
#include "Features/DynamicRupture/Typedefs.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>

#ifdef SEISSOL_KERNELS_SYCL
#include <sycl/sycl.hpp>
#endif

/**
 * What a friction law works on while it evaluates the faces of a layer, once for the host and the
 * device. The code of a law handles the points of a face in a range (see common::NumPoints): on the
 * host all points of the face at once, on the device the one point of a thread. What it keeps of
 * the points meanwhile is accordingly an array over the points of the face on the host and a single
 * value on the device (PointValues, FaultStresses<Cfg, Executor>); point() reads either at an index
 * of the range.
 */
namespace seissol::dr::friction_law {

/// What a point of the range keeps of a quantity: on the host, one value per point of the face; on
/// the device, the value of the point of the thread.
template <typename Cfg, common::RangeType Type>
using PointValues = std::conditional_t<Type == common::RangeType::CPU,
                                       std::array<Real<Cfg>, misc::NumPaddedPoints<Cfg>>,
                                       Real<Cfg>>;

/// The alignment of the values of the points of a range: that of the vectors on the host, where
/// the points of a face are processed together, that of a single value on the device.
template <typename Cfg, common::RangeType Type>
constexpr std::size_t pointAlignment() {
  if constexpr (Type == common::RangeType::CPU) {
    return Alignment;
  } else {
    return alignof(Real<Cfg>);
  }
}

/// The value of the point `index` of the range in `values`: the one of the point of the thread on
/// the device, where a single value is kept per point.
template <typename T>
SEISSOL_HOSTDEVICE constexpr T& point(T& value, std::uint32_t /*index*/) {
  return value;
}

template <typename T, std::size_t N>
SEISSOL_HOSTDEVICE constexpr T& point(std::array<T, N>& values, std::uint32_t index) {
  return values[index];
}

template <typename T, std::size_t N>
SEISSOL_HOSTDEVICE constexpr const T& point(const std::array<T, N>& values, std::uint32_t index) {
  return values[index];
}

template <typename T, std::size_t N>
SEISSOL_HOSTDEVICE constexpr T& point(T (&values)[N], std::uint32_t index) {
  return values[index];
}

template <typename T, std::size_t N>
SEISSOL_HOSTDEVICE constexpr const T& point(const T (&values)[N], std::uint32_t index) {
  return values[index];
}

/// What a friction law reads of a time step, besides the faces of the layer.
template <typename Cfg>
struct FrictionLawArgs {
  using real = Real<Cfg>;

  const FrictionLawData<Cfg>* __restrict data{nullptr};
  const real* __restrict spaceWeights{nullptr};
  const real* __restrict resampleMatrix{nullptr};
  const real* __restrict tpInverseFourierCoefficients{nullptr};
  const real* __restrict tpGridPoints{nullptr};
  const real* __restrict heatSource{nullptr};

  real timeWeights[misc::TimeSteps<Cfg>]{};
  real deltaT[misc::TimeSteps<Cfg>]{};
  real fullUpdateTime{};

  // if not null, the time at the start of the step, read on the device instead of fullUpdateTime
  const double* __restrict clock{nullptr};
};

template <typename Cfg>
struct InitialVariables {
  using real = Real<Cfg>;

  real absoluteShearTraction{};
  real localSlipRate{};
  real normalStress{};
  /// the same, before the slip rate dependent part and the clamp; the Newton solve needs it to
  /// evaluate sigma(V) itself
  real normalStressStick{};
  real stateVarReference{};
  real etaNormal{};
  /// unit slip direction; equals the normalized trial traction unless the impedance is anisotropic
  real slipDirection1{};
  real slipDirection2{};
};

/// the affine law a step of the thermal pressurization reduces to; see prepareFluidPressure
template <typename Cfg>
struct FluidPressureCoefficients {
  using real = Real<Cfg>;

  real temperatureOffset{};
  real temperatureSlope{};
  real pressureOffset{};
  real pressureSlope{};
};

/// What the rate-and-state laws and the thermal pressurization keep of a point during a time step;
/// only on the device so far, the only place where they run on the code of both
/// (FrictionLawDriver).
template <typename Cfg, common::RangeType Type>
struct LawPointState {};

template <typename Cfg>
struct LawPointState<Cfg, common::RangeType::GPU> {
  InitialVariables<Cfg> initialVariables{};
  FluidPressureCoefficients<Cfg> fluidPressure{};
};

/**
 * A friction law at the points of the range of a face: the face, the first point of the range,
 * what the law reads (data, args), and what it keeps of the points during a time step.
 */
template <typename Cfg, common::RangeType Type>
struct FrictionLawContext : public LawPointState<Cfg, Type> {
  using real = Real<Cfg>;
  static constexpr Executor Exec = common::RangeExecutor<Type>::Exec;

  std::size_t ltsFace{0};
  /// the first point of the range: the point of the thread on the device, 0 on the host
  std::uint32_t pointIndex{0};
  const FrictionLawData<Cfg>* __restrict data{nullptr};
  const FrictionLawArgs<Cfg>* __restrict args{nullptr};

  /// the time at the start of the step
  real fullUpdateTime{};

  real* __restrict sharedMemory{nullptr};
  void* item{nullptr};

  FaultStresses<Cfg, Exec> faultStresses{};
  FaultStresses<Cfg, Exec> initialStress{};
  TractionResults<Cfg, Exec> tractionResults{};
  alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> stateVariableBuffer{};
  alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> strengthBuffer{};
  /// d(strength)/d(-sigma_eff); only used for the anisotropic normal/shear coupling
  alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> strengthSlopeBuffer{};
};

/// The first point of the range of `ctx`: the point of the thread on the device; on the host, where
/// the range covers the face, 0 -- as a constant, so that the loops over the points of a face are
/// known to start at its first one.
template <typename Cfg, common::RangeType Type>
SEISSOL_HOSTDEVICE constexpr std::uint32_t firstPoint(const FrictionLawContext<Cfg, Type>& ctx) {
  if constexpr (Type == common::RangeType::CPU) {
    return 0;
  } else {
    return ctx.pointIndex;
  }
}

#ifdef __CUDACC__
template <typename Cfg>
SEISSOL_DEVICE inline void
    deviceBarrier(FrictionLawContext<Cfg, common::RangeType::GPU>& __restrict ctx) {
  __syncthreads();
}
template <typename Cfg>
SEISSOL_DEVICE inline void
    deviceWarpBarrier(FrictionLawContext<Cfg, common::RangeType::GPU>& __restrict ctx) {
  __syncwarp();
}
template <typename Cfg>
SEISSOL_DEVICE inline bool
    deviceWarpAll(FrictionLawContext<Cfg, common::RangeType::GPU>& __restrict ctx, bool value) {
  return __all_sync(0xffffffffU, static_cast<int>(value)) != 0;
}
#elif defined(__HIP__)
template <typename Cfg>
SEISSOL_DEVICE inline void
    deviceBarrier(FrictionLawContext<Cfg, common::RangeType::GPU>& __restrict ctx) {
  __syncthreads();
}
template <typename Cfg>
SEISSOL_DEVICE inline void
    deviceWarpBarrier(FrictionLawContext<Cfg, common::RangeType::GPU>& __restrict ctx) {
  // __syncwarp has no effect on current AMD GPUs (early 2026)
  // (nor does the HIP in our current CI support it)
}
template <typename Cfg>
SEISSOL_DEVICE inline bool
    deviceWarpAll(FrictionLawContext<Cfg, common::RangeType::GPU>& __restrict ctx, bool value) {
  return __all(static_cast<int>(value)) != 0;
}
#elif defined(SEISSOL_KERNELS_SYCL)
template <typename Cfg>
inline void deviceBarrier(FrictionLawContext<Cfg, common::RangeType::GPU>& __restrict ctx) {
  reinterpret_cast<sycl::nd_item<1>*>(ctx.item)->barrier(sycl::access::fence_space::local_space);
}
template <typename Cfg>
inline void deviceWarpBarrier(FrictionLawContext<Cfg, common::RangeType::GPU>& __restrict ctx) {
  const auto subgroup = reinterpret_cast<sycl::nd_item<1>*>(ctx.item)->get_sub_group();
  sycl::group_barrier(subgroup);
}
template <typename Cfg>
inline bool deviceWarpAll(FrictionLawContext<Cfg, common::RangeType::GPU>& __restrict ctx,
                          bool value) {
  const auto subgroup = reinterpret_cast<sycl::nd_item<1>*>(ctx.item)->get_sub_group();
  return sycl::all_of_group(subgroup, value);
}
#else
template <typename Cfg>
inline void deviceBarrier(FrictionLawContext<Cfg, common::RangeType::GPU>& __restrict /*ctx*/) {}
template <typename Cfg>
inline void deviceWarpBarrier(FrictionLawContext<Cfg, common::RangeType::GPU>& __restrict /*ctx*/) {
}
template <typename Cfg>
inline bool deviceWarpAll(FrictionLawContext<Cfg, common::RangeType::GPU>& __restrict /*ctx*/,
                          bool /*value*/) {
  return true;
}
#endif

/// The value `toResample` of the point of the thread, resampled together with the values of the
/// other points of its face (see resample); on the device.
template <typename Cfg>
SEISSOL_DEVICE inline Real<Cfg>
    resampleVariable(FrictionLawContext<Cfg, common::RangeType::GPU>& __restrict ctx,
                     Real<Cfg> toResample) {
  constexpr auto Dim0 = misc::dimSize<init::resample<Cfg>, 0>();
  constexpr auto Dim1 = misc::dimSize<init::resample<Cfg>, 1>();
  static_assert(Dim0 == misc::NumPaddedPointsSingleSim<Cfg>);
  static_assert(Dim0 >= Dim1);

  ctx.sharedMemory[ctx.pointIndex] = toResample;
  deviceBarrier(ctx);

  const auto simPointIndex = ctx.pointIndex / Cfg::NumSimulations;
  const auto simId = ctx.pointIndex % Cfg::NumSimulations;
  constexpr uint32_t SimPointStride =
      multisim::MultisimHelperWrapper<Cfg>::MultisimEnabled ? Dim1 : 1U;
  constexpr uint32_t DataPointStride =
      multisim::MultisimHelperWrapper<Cfg>::MultisimEnabled ? 1U : Dim0;

  Real<Cfg> result{0};
  for (uint32_t i = 0; i < Dim1; ++i) {
    result += ctx.args->resampleMatrix[simPointIndex * SimPointStride + i * DataPointStride] *
              ctx.sharedMemory[i * Cfg::NumSimulations + simId];
  }
  deviceBarrier(ctx);

  return result;
}

/// The values `values` at the points of a face resampled with the matrix `resampleMatrix`; on the
/// host.
template <typename Cfg>
inline void resampleOnHost(const Real<Cfg>* resampleMatrix,
                           const Real<Cfg> (&values)[misc::NumPaddedPoints<Cfg>],
                           std::array<Real<Cfg>, misc::NumPaddedPoints<Cfg>>& resampled) {
  dynamicRupture::kernel::resampleParameter<Cfg> resampleKrnl;
  resampleKrnl.resample = resampleMatrix;
  resampleKrnl.originalQ = values;
  resampleKrnl.resampledQ = resampled.data();
  resampleKrnl.execute();
}

/**
 * The values `values` at the points of the face of `ctx`, projected onto the polynomials of the
 * face of degree less than the order of the configuration and evaluated at its points again; the
 * ones of the points of the range go to `resampled`. All points of the face resample together: on
 * the device, the threads of a face meet in its shared memory.
 */
template <typename Cfg, common::RangeType Type>
SEISSOL_HOSTDEVICE inline void resample(FrictionLawContext<Cfg, Type>& ctx,
                                        const Real<Cfg> (&values)[misc::NumPaddedPoints<Cfg>],
                                        PointValues<Cfg, Type>& resampled) {
  if constexpr (Type == common::RangeType::CPU) {
    resampleOnHost<Cfg>(ctx.args->resampleMatrix, values, resampled);
  } else {
    resampled = resampleVariable(ctx, values[ctx.pointIndex]);
  }
}

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_FRICTIONLAWCONTEXT_H_
