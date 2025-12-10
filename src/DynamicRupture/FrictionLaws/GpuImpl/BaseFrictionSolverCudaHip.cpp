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

template <typename Cfg>
constexpr std::size_t PaddedMultiple =
    safeblockMultiple(seissol::dr::misc::NumPaddedPoints<Cfg>, BlockTargetsize);

template <typename Cfg, typename T>
__launch_bounds__(PaddedMultiple<Cfg>* seissol::dr::misc::NumPaddedPoints<Cfg>) __global__
    void flkernelwrapper(std::size_t elements, FrictionLawArgs<Cfg> args) {
  using real = Real<Cfg>;
  FrictionLawContext<Cfg> ctx{};

  ctx.data = args.data;
  ctx.args = &args;

  __shared__ real shm[PaddedMultiple<Cfg> * seissol::dr::misc::NumPaddedPoints<Cfg>];
  ctx.sharedMemory = &shm[threadIdx.z * seissol::dr::misc::NumPaddedPoints<Cfg>];
  // ctx.item = nullptr;

  ctx.ltsFace = blockIdx.x * PaddedMultiple<Cfg> + threadIdx.z;
  ctx.pointIndex = threadIdx.x + threadIdx.y * multisim::NumSimulations<Cfg>;

  if (ctx.ltsFace < elements) {
    seissol::dr::friction_law::gpu::BaseFrictionSolver<Cfg, T>::evaluatePoint(ctx);
  }
}
} // namespace

template <typename Cfg, typename T>
void BaseFrictionSolver<Cfg, T>::evaluateKernel(seissol::parallel::runtime::StreamRuntime& runtime,
                                                double fullUpdateTime,
                                                const double* timeWeights,
                                                const FrictionSolver::FrictionTime& frictionTime) {
#ifdef __CUDACC__
  using StreamT = cudaStream_t;
#endif
#ifdef __HIP__
  using StreamT = hipStream_t;
#endif
  auto stream = reinterpret_cast<StreamT>(runtime.stream());
  dim3 block(
      multisim::NumSimulations<Cfg>, misc::NumPaddedPointsSingleSim<Cfg>, PaddedMultiple<Cfg>);
  dim3 grid((this->currLayerSize + PaddedMultiple<Cfg> - 1) / PaddedMultiple<Cfg>);

  FrictionLawArgs<Cfg> args{};
  args.data = this->data;
  args.spaceWeights = this->devSpaceWeights;
  args.resampleMatrix = this->resampleMatrix;
  args.tpInverseFourierCoefficients = this->devTpInverseFourierCoefficients;
  args.tpGridPoints = this->devTpGridPoints;
  args.heatSource = this->devHeatSource;
  std::copy_n(timeWeights, misc::TimeSteps<Cfg>, args.timeWeights);
  std::copy_n(frictionTime.deltaT.data(), misc::TimeSteps<Cfg>, args.deltaT);
  args.sumDt = frictionTime.sumDt;
  args.fullUpdateTime = fullUpdateTime;

  flkernelwrapper<Cfg, T><<<grid, block, 0, stream>>>(this->currLayerSize, args);
}

#define SEISSOL_CONFIGITER(cfg)                                                                    \
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
