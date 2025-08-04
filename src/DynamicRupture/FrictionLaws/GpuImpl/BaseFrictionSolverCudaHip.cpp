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
#include <Common/Constants.h>

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
    safeblockMultiple(seissol::dr::misc::NumPaddedPoints<Cfg>, BlockTargetsize);

template <typename T>
__launch_bounds__(PaddedMultiple* seissol::dr::misc::NumPaddedPoints<Cfg>) __global__
    void flkernelwrapper(std::size_t elements, FrictionLawArgs args) {
  FrictionLawContext ctx{};

  ctx.data = args.data;
  ctx.args = &args;

  __shared__ real shm[PaddedMultiple * seissol::dr::misc::NumPaddedPoints<Cfg>];
  ctx.sharedMemory = &shm[threadIdx.z * seissol::dr::misc::NumPaddedPoints<Cfg>];
  // ctx.item = nullptr;

  ctx.ltsFace = blockIdx.x * PaddedMultiple + threadIdx.z;
  ctx.pointIndex = threadIdx.x + threadIdx.y * multisim::NumSimulations;

  if (ctx.ltsFace < elements) {
    seissol::dr::friction_law::gpu::BaseFrictionSolver<T>::evaluatePoint(ctx);
  }
}
} // namespace

namespace seissol::dr::friction_law::gpu {

template <typename T>
void BaseFrictionSolver<T>::evaluateKernel(seissol::parallel::runtime::StreamRuntime& runtime,
                                           double fullUpdateTime,
                                           const double* timeWeights,
                                           const FrictionTime& frictionTime) {
#ifdef __CUDACC__
  using StreamT = cudaStream_t;
#endif
#ifdef __HIP__
  using StreamT = hipStream_t;
#endif
  auto stream = reinterpret_cast<StreamT>(runtime.stream());
  dim3 block(multisim::NumSimulations, misc::NumPaddedPointsSingleSim<Cfg>, PaddedMultiple);
  dim3 grid((this->currLayerSize + PaddedMultiple - 1) / PaddedMultiple);

  FrictionLawArgs args{};
  args.data = data;
  args.spaceWeights = devSpaceWeights;
  args.resampleMatrix = resampleMatrix;
  args.tpInverseFourierCoefficients = devTpInverseFourierCoefficients;
  args.tpGridPoints = devTpGridPoints;
  args.heatSource = devHeatSource;
  std::copy_n(timeWeights, Cfg::ConvergenceOrder, args.timeWeights);
  std::copy_n(frictionTime.deltaT.data(), Cfg::ConvergenceOrder, args.deltaT);
  args.sumDt = frictionTime.sumDt;
  args.fullUpdateTime = fullUpdateTime;

  flkernelwrapper<T><<<grid, block, 0, stream>>>(this->currLayerSize, args);
}

#define _H_(cfg)                                                                                   \
  template class BaseFrictionSolver<cfg, NoFault<cfg>>;                                            \
  template class BaseFrictionSolver<                                                               \
      cfg,                                                                                         \
      LinearSlipWeakeningBase<cfg, LinearSlipWeakeningLaw<cfg, NoSpecialization<cfg>>>>;           \
  template class BaseFrictionSolver<                                                               \
      cfg,                                                                                         \
      LinearSlipWeakeningBase<cfg, LinearSlipWeakeningLaw<cfg, BiMaterialFault<cfg>>>>;            \
  template class BaseFrictionSolver<                                                               \
      cfg,                                                                                         \
      LinearSlipWeakeningBase<cfg, LinearSlipWeakeningLaw<cfg, TPApprox<cfg>>>>;                   \
  template class BaseFrictionSolver<                                                               \
      cfg,                                                                                         \
      RateAndStateBase<cfg,                                                                        \
                       SlowVelocityWeakeningLaw<cfg, AgingLaw<cfg, NoTP<cfg>>, NoTP<cfg>>,         \
                       NoTP<cfg>>>;                                                                \
  template class BaseFrictionSolver<                                                               \
      cfg,                                                                                         \
      RateAndStateBase<cfg,                                                                        \
                       SlowVelocityWeakeningLaw<cfg, SlipLaw<cfg, NoTP<cfg>>, NoTP<cfg>>,          \
                       NoTP<cfg>>>;                                                                \
  template class BaseFrictionSolver<                                                               \
      cfg,                                                                                         \
      RateAndStateBase<cfg, FastVelocityWeakeningLaw<cfg, NoTP<cfg>>, NoTP<cfg>>>;                 \
  template class BaseFrictionSolver<                                                               \
      cfg,                                                                                         \
      RateAndStateBase<cfg, SevereVelocityWeakeningLaw<cfg, NoTP<cfg>>, NoTP<cfg>>>;               \
  template class BaseFrictionSolver<                                                               \
      cfg,                                                                                         \
      RateAndStateBase<cfg,                                                                        \
                       SlowVelocityWeakeningLaw<cfg,                                               \
                                                AgingLaw<cfg, ThermalPressurization<cfg>>,         \
                                                ThermalPressurization<cfg>>,                       \
                       ThermalPressurization<cfg>>>;                                               \
  template class BaseFrictionSolver<                                                               \
      cfg,                                                                                         \
      RateAndStateBase<cfg,                                                                        \
                       SlowVelocityWeakeningLaw<cfg,                                               \
                                                SlipLaw<cfg, ThermalPressurization<cfg>>,          \
                                                ThermalPressurization<cfg>>,                       \
                       ThermalPressurization<cfg>>>;                                               \
  template class BaseFrictionSolver<                                                               \
      cfg,                                                                                         \
      RateAndStateBase<cfg,                                                                        \
                       FastVelocityWeakeningLaw<cfg, ThermalPressurization<cfg>>,                  \
                       ThermalPressurization<cfg>>>;                                               \
  template class BaseFrictionSolver<                                                               \
      cfg,                                                                                         \
      RateAndStateBase<cfg,                                                                        \
                       SevereVelocityWeakeningLaw<cfg, ThermalPressurization<cfg>>,                \
                       ThermalPressurization<cfg>>>;                                               \
  template class BaseFrictionSolver<cfg, ImposedSlipRates<cfg, YoffeSTF<cfg>>>;                    \
  template class BaseFrictionSolver<cfg, ImposedSlipRates<cfg, GaussianSTF<cfg>>>;                 \
  template class BaseFrictionSolver<cfg, ImposedSlipRates<cfg, DeltaSTF<cfg>>>;

#include "ConfigInclude.h"

} // namespace seissol::dr::friction_law::gpu
