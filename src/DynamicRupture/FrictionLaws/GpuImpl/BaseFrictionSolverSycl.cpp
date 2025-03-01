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

#include <sycl/sycl.hpp>

namespace seissol::dr::friction_law::gpu {

template <typename T>
void BaseFrictionSolver<T>::evaluateKernel(seissol::parallel::runtime::StreamRuntime& runtime,
                                           real fullUpdateTime) {
  auto* data{this->data};
  auto* devTimeWeights{this->devTimeWeights};
  auto* devSpaceWeights{this->devSpaceWeights};
  auto* resampleMatrix{this->resampleMatrix};
  auto devFullUpdateTime{fullUpdateTime};

  auto* TpInverseFourierCoefficients{this->devTpInverseFourierCoefficients};
  auto* TpGridPoints{this->devTpGridPoints};
  auto* HeatSource{this->devHeatSource};

  auto* queue = reinterpret_cast<sycl::queue*>(runtime.stream());

  sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints}, {misc::NumPaddedPoints}};
  queue->submit([&](sycl::handler& cgh) {
    // NOLINTNEXTLINE
    sycl::accessor<real, 1, sycl::access::mode::read_write, sycl::access::target::local>
        sharedMemory(misc::NumPaddedPoints, cgh);

    cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
      FrictionLawContext ctx{};
      ctx.sharedMemory = &sharedMemory[0];
      ctx.item = reinterpret_cast<void*>(&item);
      ctx.data = data;
      ctx.devTimeWeights = devTimeWeights;
      ctx.devSpaceWeights = devSpaceWeights;
      ctx.resampleMatrix = resampleMatrix;
      ctx.fullUpdateTime = devFullUpdateTime;
      ctx.TpInverseFourierCoefficients = TpInverseFourierCoefficients;
      ctx.TpGridPoints = TpGridPoints;
      ctx.HeatSource = HeatSource;

      const auto ltsFace = item.get_group().get_group_id(0);
      const auto pointIndex = item.get_local_id(0);

      ctx.ltsFace = ltsFace;
      ctx.pointIndex = pointIndex;

      evaluatePoint(ctx);
    });
  });
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
