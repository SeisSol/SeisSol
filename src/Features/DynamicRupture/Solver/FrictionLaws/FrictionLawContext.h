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

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>

#ifdef SEISSOL_KERNELS_SYCL
#include <sycl/sycl.hpp>
#endif

/**
 * What a friction law works on while it evaluates the faces of a layer, once for the host and the
 * device. The code of a law handles the points of a face in a range (forEachPoint): on the host all
 * points of the face at once, on the device the one point of a thread. What it keeps of the points
 * meanwhile is accordingly an array over the points of the face on the host and the value of the
 * point of the thread on the device (PointValues, FaultStresses<Cfg, Executor>); either is indexed
 * with the index of the point in the face, as the storage of the face is.
 */
namespace seissol::dr::friction_law {

/// What the range keeps of a quantity of the type T at its points: on the host, one value per point
/// of the face; on the device, the value of the point of the thread. Indexed with the index of the
/// point in the face.
template <typename T, typename Cfg, common::RangeType Type>
using PointArray = std::conditional_t<Type == common::RangeType::CPU,
                                      std::array<T, misc::NumPaddedPoints<Cfg>>,
                                      ThreadPointValue<T>>;

/// What the range keeps of a real quantity at its points (see PointArray).
template <typename Cfg, common::RangeType Type>
using PointValues = PointArray<Real<Cfg>, Cfg, Type>;

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

/**
 * A friction law at the points of the range of a face: the face, the first point of the range,
 * what the law reads (data, args), and what it keeps of the points during a time step.
 */
template <typename Cfg, common::RangeType Type>
struct FrictionLawContext {
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

/**
 * Runs `body(pointIndex)` for the points of the range of `ctx`, with the index of a point in its
 * face: on the host for all points of the face, as one loop that is vectorized unless `Vectorized`
 * is false (the points are independent of each other); on the device for the point of the thread.
 * The index addresses both the storage of the face (ctx.data->...[ctx.ltsFace][pointIndex]) and
 * what the range keeps of its points (PointValues, ctx.faultStresses, ...).
 */
template <bool Vectorized = true, typename Cfg, common::RangeType Type, typename BodyT>
SEISSOL_HOSTDEVICE inline void forEachPoint(const FrictionLawContext<Cfg, Type>& ctx,
                                            const BodyT& body) {
  using Range = common::NumPoints<Cfg, Type>::Range;
  if constexpr (Vectorized) {
#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (auto index = Range::Start; index < Range::End; index += Range::Step) {
      body(firstPoint(ctx) + index);
    }
  } else {
    for (auto index = Range::Start; index < Range::End; index += Range::Step) {
      body(firstPoint(ctx) + index);
    }
  }
}

/**
 * forEachPoint without vectorizing the loop over the points on the host: for code that branches by
 * the point and calls transcendental functions, all of whose branches every point would evaluate in
 * a vector (e.g. the source time functions of ImposedSlipRates).
 */
template <typename Cfg, common::RangeType Type, typename BodyT>
SEISSOL_HOSTDEVICE inline void forEachPointScalar(const FrictionLawContext<Cfg, Type>& ctx,
                                                  const BodyT& body) {
  forEachPoint<false>(ctx, body);
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

/**
 * Whether `flags` holds at all points of the range of `ctx`: on the host at all points of the face;
 * on the device at the points of the warp of the thread, whose threads thus go on or stop together.
 */
template <typename Cfg, common::RangeType Type>
SEISSOL_HOSTDEVICE inline bool allPoints(FrictionLawContext<Cfg, Type>& ctx,
                                         const PointArray<std::int32_t, Cfg, Type>& flags) {
  if constexpr (Type == common::RangeType::CPU) {
    return std::all_of(flags.begin(), flags.end(), [](std::int32_t flag) { return flag != 0; });
  } else {
    deviceWarpBarrier(ctx);
    return deviceWarpAll(ctx, flags[ctx.pointIndex] != 0);
  }
}

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
                           const Real<Cfg>* values,
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
    resampled[ctx.pointIndex] = resampleVariable(ctx, values[ctx.pointIndex]);
  }
}

/// resample for the values of the points of the range (see PointValues)
template <typename Cfg, common::RangeType Type>
SEISSOL_HOSTDEVICE inline void resample(FrictionLawContext<Cfg, Type>& ctx,
                                        const PointValues<Cfg, Type>& values,
                                        PointValues<Cfg, Type>& resampled) {
  if constexpr (Type == common::RangeType::CPU) {
    resampleOnHost<Cfg>(ctx.args->resampleMatrix, values.data(), resampled);
  } else {
    resampled[ctx.pointIndex] = resampleVariable(ctx, values[ctx.pointIndex]);
  }
}

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_FRICTIONLAWCONTEXT_H_
