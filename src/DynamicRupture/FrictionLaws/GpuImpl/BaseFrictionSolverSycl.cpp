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
#include "ThermalPressurization<Cfg>/NoTP<Cfg>.h"
#include "ThermalPressurization<Cfg>/ThermalPressurization<Cfg>.h"

#include <sycl/sycl.hpp>

namespace seissol::dr::friction_law::gpu {

template <typename T>
void BaseFrictionSolver<Cfg, T>::evaluateKernel(seissol::parallel::runtime::StreamRuntime& runtime,
                                                double fullUpdateTime,
                                                const double* timeWeights,
                                                const FrictionTime& frictionTime) {
  auto* queue = reinterpret_cast<sycl::queue*>(runtime.stream());

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

  sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints<Cfg>},
                     {misc::NumPaddedPoints<Cfg>}};
  queue->submit([&](sycl::handler& cgh) {
    // NOLINTNEXTLINE
    sycl::local_accessor<real> sharedMemory(misc::NumPaddedPoints<Cfg>, cgh);

    cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
      FrictionLawContext ctx{};
      ctx.sharedMemory = &sharedMemory[0];
      ctx.item = reinterpret_cast<void*>(&item);
      ctx.data = args.data;
      ctx.args = &args;

      const auto ltsFace = item.get_group().get_group_id(0);
      const auto pointIndex = item.get_local_id(0);

      ctx.ltsFace = ltsFace;
      ctx.pointIndex = pointIndex;

      evaluatePoint(ctx);
    });
  });
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
