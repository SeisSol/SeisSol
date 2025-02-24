// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "BaseFrictionSolver.h"

#include "AgingLaw.h"
#include "FastVelocityWeakeningLaw.h"
#include "FrictionSolverInterface.h"
#include "ImposedSlipRates.h"
#include "LinearSlipWeakening.h"
#include "NoFault.h"
#include "RateAndState.h"
#include "SevereVelocityWeakeningLaw.h"
#include "SlipLaw.h"
#include "SlowVelocityWeakeningLaw.h"
#include "SourceTimeFunction.h"
#include "ThermalPressurization/NoTP.h"
#include "ThermalPressurization/ThermalPressurization.h"

#ifdef __HIP__
#include "hip/hip_runtime.h"
#endif

namespace {
using namespace seissol::dr::friction_law::gpu;

constexpr std::size_t safeblockMultiple(std::size_t block, std::size_t maxmult) {
  const auto unsafe = (maxmult + block - 1) / block;
  if (unsafe == 0) {
    return 1;
  } else {
    return unsafe;
  }
}

constexpr std::size_t BlockTargetsize = 256;
constexpr std::size_t PaddedMultiple =
    safeblockMultiple(seissol::dr::misc::NumPaddedPoints, BlockTargetsize);

template <typename T>
__launch_bounds__(PaddedMultiple* seissol::dr::misc::NumPaddedPoints) __global__
    void flkernelwrapper(int elements,
                         FrictionLawData* data,
                         double* timeWeights,
                         real* spaceWeights,
                         real* resampleMatrix,
                         real fullUpdateTime,
                         real* TpInverseFourierCoefficients,
                         real* TpGridPoints,
                         real* HeatSource) {
  FrictionLawContext ctx{};

  ctx.data = data;
  ctx.devTimeWeights = timeWeights;
  ctx.devSpaceWeights = spaceWeights;
  ctx.resampleMatrix = resampleMatrix;
  ctx.fullUpdateTime = fullUpdateTime;

  ctx.TpInverseFourierCoefficients = TpInverseFourierCoefficients;
  ctx.TpGridPoints = TpGridPoints;
  ctx.HeatSource = HeatSource;

  __shared__ real shm[PaddedMultiple * seissol::dr::misc::NumPaddedPoints];
  ctx.sharedMemory = &shm[threadIdx.y * seissol::dr::misc::NumPaddedPoints];
  // ctx.item = nullptr;

  ctx.ltsFace = blockIdx.x * blockDim.y + threadIdx.y;
  ctx.pointIndex = threadIdx.x;

  if (ctx.ltsFace < elements) {
    seissol::dr::friction_law::gpu::BaseFrictionSolver<T>::evaluatePoint(ctx);
  }
}
} // namespace

namespace seissol::dr::friction_law::gpu {

template <typename T>
void BaseFrictionSolver<T>::evaluateKernel(seissol::parallel::runtime::StreamRuntime& runtime,
                                           real fullUpdateTime) {
#ifdef __CUDACC__
  using StreamT = cudaStream_t;
#endif
#ifdef __HIP__
  using StreamT = hipStream_t;
#endif
  auto stream = reinterpret_cast<StreamT>(runtime.stream());
  dim3 block(misc::NumPaddedPoints, PaddedMultiple);
  dim3 grid((this->currLayerSize + PaddedMultiple - 1) / PaddedMultiple);
  flkernelwrapper<T><<<grid, block, 0, stream>>>(this->currLayerSize,
                                                 data,
                                                 devTimeWeights,
                                                 devSpaceWeights,
                                                 resampleMatrix,
                                                 fullUpdateTime,
                                                 devTpInverseFourierCoefficients,
                                                 devTpGridPoints,
                                                 devHeatSource);
}

template class BaseFrictionSolver<NoFault>;
template class BaseFrictionSolver<
    LinearSlipWeakeningBase<LinearSlipWeakeningLaw<NoSpecialization>>>;
template class BaseFrictionSolver<LinearSlipWeakeningBase<LinearSlipWeakeningLaw<BiMaterialFault>>>;
template class BaseFrictionSolver<LinearSlipWeakeningBase<LinearSlipWeakeningLaw<TPApprox>>>;
template class BaseFrictionSolver<
    RateAndStateBase<SlowVelocityWeakeningLaw<AgingLaw<NoTP>, NoTP>, NoTP>>;
template class BaseFrictionSolver<
    RateAndStateBase<SlowVelocityWeakeningLaw<SlipLaw<NoTP>, NoTP>, NoTP>>;
template class BaseFrictionSolver<RateAndStateBase<FastVelocityWeakeningLaw<NoTP>, NoTP>>;
template class BaseFrictionSolver<RateAndStateBase<SevereVelocityWeakeningLaw<NoTP>, NoTP>>;
template class BaseFrictionSolver<RateAndStateBase<
    SlowVelocityWeakeningLaw<AgingLaw<ThermalPressurization>, ThermalPressurization>,
    ThermalPressurization>>;
template class BaseFrictionSolver<RateAndStateBase<
    SlowVelocityWeakeningLaw<SlipLaw<ThermalPressurization>, ThermalPressurization>,
    ThermalPressurization>>;
template class BaseFrictionSolver<
    RateAndStateBase<FastVelocityWeakeningLaw<ThermalPressurization>, ThermalPressurization>>;
template class BaseFrictionSolver<
    RateAndStateBase<SevereVelocityWeakeningLaw<ThermalPressurization>, ThermalPressurization>>;
template class BaseFrictionSolver<ImposedSlipRates<YoffeSTF>>;
template class BaseFrictionSolver<ImposedSlipRates<GaussianSTF>>;
template class BaseFrictionSolver<ImposedSlipRates<DeltaSTF>>;

} // namespace seissol::dr::friction_law::gpu
