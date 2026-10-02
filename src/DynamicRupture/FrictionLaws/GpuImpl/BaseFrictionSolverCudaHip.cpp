// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "AgingLaw.h"
#include "BaseFrictionSolver.h"
#include "Common/Real.h"
#include "Config.h"
#include "DynamicRupture/Misc.h"
#include "FastVelocityWeakeningLaw.h"
#include "ImposedSlipRates.h"
#include "LinearSlipWeakening.h"
#include "NoFault.h"
#include "Parallel/Runtime/Stream.h"
#include "RateAndState.h"
#include "SevereVelocityWeakeningLaw.h"
#include "SlipLaw.h"
#include "SlowVelocityWeakeningLaw.h"
#include "Solver/MultipleSimulations.h"
#include "SourceTimeFunction.h"
#include "ThermalPressurization/NoTP.h"
#include "ThermalPressurization/ThermalPressurization.h"

#include <cstddef>

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

// __grid_constant__ needs sm_70 or higher
#if defined(__CUDACC__) && defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 700
#define SEISSOL_GRID_CONSTANT __grid_constant__
#else
#define SEISSOL_GRID_CONSTANT
#endif

template <typename Cfg, typename T>
__launch_bounds__(PaddedMultiple<Cfg>* seissol::dr::misc::NumPaddedPoints<Cfg>) __global__
    void flkernelwrapper(const std::size_t elements,
                         const SEISSOL_GRID_CONSTANT FrictionLawArgs<Cfg> args) {
  FrictionLawContext<Cfg> ctx{};

  ctx.data = args.data;
  ctx.args = &args;

  __shared__ Real<Cfg> shm[PaddedMultiple<Cfg> * seissol::dr::misc::NumPaddedPoints<Cfg>];
  ctx.sharedMemory = &shm[threadIdx.z * seissol::dr::misc::NumPaddedPoints<Cfg>];
  // ctx.item = nullptr;

  ctx.ltsFace = blockIdx.x * PaddedMultiple<Cfg> + threadIdx.z;
  ctx.pointIndex = threadIdx.x + threadIdx.y * Cfg::NumSimulations;

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
  dim3 block(Cfg::NumSimulations, misc::NumPaddedPointsSingleSim<Cfg>, PaddedMultiple<Cfg>);
  dim3 grid((this->currLayerSize_ + PaddedMultiple<Cfg> - 1) / PaddedMultiple<Cfg>);

  FrictionLawArgs<Cfg> args{};
  args.data = this->data_;
  args.spaceWeights = this->devSpaceWeights_;
  args.resampleMatrix = this->resampleMatrix_;
  args.tpInverseFourierCoefficients = this->devTpInverseFourierCoefficients_;
  args.tpGridPoints = this->devTpGridPoints_;
  args.heatSource = this->devHeatSource_;
  std::copy_n(timeWeights, misc::TimeSteps<Cfg>, args.timeWeights);
  std::copy_n(frictionTime.deltaT.data(), misc::TimeSteps<Cfg>, args.deltaT);
  args.fullUpdateTime = fullUpdateTime;

  flkernelwrapper<Cfg, T><<<grid, block, 0, stream>>>(this->currLayerSize_, args);
}

#define SEISSOL_INSTANTIATE(Cfg)                                                                   \
  template class BaseFrictionSolver<Cfg, NoFault<Cfg>>;                                            \
  template class BaseFrictionSolver<                                                               \
      Cfg,                                                                                         \
      LinearSlipWeakeningBase<Cfg, LinearSlipWeakeningLaw<Cfg, NoSpecialization<Cfg>>>>;           \
  template class BaseFrictionSolver<                                                               \
      Cfg,                                                                                         \
      LinearSlipWeakeningBase<Cfg, LinearSlipWeakeningLaw<Cfg, BiMaterialFault<Cfg>>>>;            \
  template class BaseFrictionSolver<                                                               \
      Cfg,                                                                                         \
      LinearSlipWeakeningBase<Cfg, LinearSlipWeakeningLaw<Cfg, TPApprox<Cfg>>>>;                   \
  template class BaseFrictionSolver<                                                               \
      Cfg,                                                                                         \
      RateAndStateBase<Cfg,                                                                        \
                       SlowVelocityWeakeningLaw<Cfg, AgingLaw<Cfg, NoTP<Cfg>>, NoTP<Cfg>>,         \
                       NoTP<Cfg>>>;                                                                \
  template class BaseFrictionSolver<                                                               \
      Cfg,                                                                                         \
      RateAndStateBase<Cfg,                                                                        \
                       SlowVelocityWeakeningLaw<Cfg, SlipLaw<Cfg, NoTP<Cfg>>, NoTP<Cfg>>,          \
                       NoTP<Cfg>>>;                                                                \
  template class BaseFrictionSolver<                                                               \
      Cfg,                                                                                         \
      RateAndStateBase<Cfg, FastVelocityWeakeningLaw<Cfg, NoTP<Cfg>>, NoTP<Cfg>>>;                 \
  template class BaseFrictionSolver<                                                               \
      Cfg,                                                                                         \
      RateAndStateBase<Cfg, SevereVelocityWeakeningLaw<Cfg, NoTP<Cfg>>, NoTP<Cfg>>>;               \
  template class BaseFrictionSolver<                                                               \
      Cfg,                                                                                         \
      RateAndStateBase<Cfg,                                                                        \
                       SlowVelocityWeakeningLaw<Cfg,                                               \
                                                AgingLaw<Cfg, ThermalPressurization<Cfg>>,         \
                                                ThermalPressurization<Cfg>>,                       \
                       ThermalPressurization<Cfg>>>;                                               \
  template class BaseFrictionSolver<                                                               \
      Cfg,                                                                                         \
      RateAndStateBase<Cfg,                                                                        \
                       SlowVelocityWeakeningLaw<Cfg,                                               \
                                                SlipLaw<Cfg, ThermalPressurization<Cfg>>,          \
                                                ThermalPressurization<Cfg>>,                       \
                       ThermalPressurization<Cfg>>>;                                               \
  template class BaseFrictionSolver<                                                               \
      Cfg,                                                                                         \
      RateAndStateBase<Cfg,                                                                        \
                       FastVelocityWeakeningLaw<Cfg, ThermalPressurization<Cfg>>,                  \
                       ThermalPressurization<Cfg>>>;                                               \
  template class BaseFrictionSolver<                                                               \
      Cfg,                                                                                         \
      RateAndStateBase<Cfg,                                                                        \
                       SevereVelocityWeakeningLaw<Cfg, ThermalPressurization<Cfg>>,                \
                       ThermalPressurization<Cfg>>>;                                               \
  template class BaseFrictionSolver<Cfg, ImposedSlipRates<Cfg, YoffeSTF<Cfg>>>;                    \
  template class BaseFrictionSolver<Cfg, ImposedSlipRates<Cfg, GaussianSTF<Cfg>>>;                 \
  template class BaseFrictionSolver<Cfg, ImposedSlipRates<Cfg, DeltaSTF<Cfg>>>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::dr::friction_law::gpu
