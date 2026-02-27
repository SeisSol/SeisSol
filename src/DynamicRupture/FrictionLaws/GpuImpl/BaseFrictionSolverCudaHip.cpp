// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "AgingLaw.h"
#include "BaseFrictionSolver.h"
#include "Common/Constants.h"
#include "DynamicRupture/Misc.h"
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

namespace seissol::dr::friction_law::gpu {

namespace {

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
    void flkernelwrapper(std::size_t elements, FrictionLawArgs args) {
  FrictionLawContext ctx{};

  ctx.data = args.data;
  ctx.args = &args;

  __shared__ real shm[PaddedMultiple * seissol::dr::misc::NumPaddedPoints];
  ctx.sharedMemory = &shm[threadIdx.z * seissol::dr::misc::NumPaddedPoints];
  // ctx.item = nullptr;

  ctx.ltsFace = blockIdx.x * PaddedMultiple + threadIdx.z;
  ctx.pointIndex = threadIdx.x + threadIdx.y * multisim::NumSimulations;

  if (ctx.ltsFace < elements) {
    seissol::dr::friction_law::gpu::BaseFrictionSolver<T>::evaluatePoint(ctx);
  }
}
} // namespace

template <typename T>
void BaseFrictionSolver<T>::evaluateKernel(seissol::parallel::runtime::StreamRuntime& runtime,
                                           real fullUpdateTime,
                                           const double* timeWeights,
                                           const FrictionTime& frictionTime) {
#ifdef __CUDACC__
  using StreamT = cudaStream_t;
#endif
#ifdef __HIP__
  using StreamT = hipStream_t;
#endif
  auto stream = reinterpret_cast<StreamT>(runtime.stream());
  dim3 block(multisim::NumSimulations, misc::NumPaddedPointsSingleSim, PaddedMultiple);
  dim3 grid((this->currLayerSize + PaddedMultiple - 1) / PaddedMultiple);

  FrictionLawArgs args{};
  args.data = data;
  args.spaceWeights = devSpaceWeights;
  args.resampleMatrix = resampleMatrix;
  args.tpInverseFourierCoefficients = devTpInverseFourierCoefficients;
  args.tpGridPoints = devTpGridPoints;
  args.heatSource = devHeatSource;
  std::copy_n(timeWeights, misc::TimeSteps, args.timeWeights);
  std::copy_n(frictionTime.deltaT.data(), misc::TimeSteps, args.deltaT);
  args.fullUpdateTime = fullUpdateTime;

  flkernelwrapper<T><<<grid, block, 0, stream>>>(this->currLayerSize, args);
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
