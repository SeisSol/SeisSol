// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_

#include "Common/Constants.h"
#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverDetails.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/Misc.h"
#include "FrictionSolverInterface.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Numerical/Functions.h"
#include <algorithm>

#include "Common/Marker.h"

#ifdef SEISSOL_KERNELS_SYCL
#include <sycl/sycl.hpp>
#endif

namespace seissol::dr::friction_law::gpu {
struct InitialVariables {
  real absoluteShearTraction{};
  real localSlipRate{};
  real normalStress{};
  real stateVarReference{};
};

struct FrictionLawArgs {
  const FrictionLawData* __restrict data{nullptr};
  const real* __restrict spaceWeights{nullptr};
  const real* __restrict resampleMatrix{nullptr};
  const real* __restrict tpInverseFourierCoefficients{nullptr};
  const real* __restrict tpGridPoints{nullptr};
  const real* __restrict heatSource{nullptr};

  real fullUpdateTime;
  double timeWeights[misc::TimeSteps];
  real deltaT[misc::TimeSteps];
  real sumDt;
};

struct FrictionLawContext {
  std::size_t ltsFace;
  std::uint32_t pointIndex;
  const FrictionLawData* __restrict data;
  const FrictionLawArgs* __restrict args;

  real* __restrict sharedMemory;
  void* item;

  FaultStresses<Executor::Device> faultStresses{};
  TractionResults<Executor::Device> tractionResults{};
  real stateVariableBuffer{};
  real strengthBuffer{};
  InitialVariables initialVariables{};
};

#ifdef __CUDACC__
SEISSOL_DEVICE inline void deviceBarrier(FrictionLawContext& ctx) { __syncthreads(); }
#elif defined(__HIP__)
SEISSOL_DEVICE inline void deviceBarrier(FrictionLawContext& ctx) { __syncthreads(); }
#elif defined(SEISSOL_KERNELS_SYCL)
inline void deviceBarrier(FrictionLawContext& ctx) {
  reinterpret_cast<sycl::nd_item<1>*>(ctx.item)->barrier(sycl::access::fence_space::local_space);
}
#else
inline void deviceBarrier(FrictionLawContext& ctx) {}
#endif

template <typename Derived>
class BaseFrictionSolver : public FrictionSolverDetails {
  public:
  explicit BaseFrictionSolver<Derived>(seissol::initializer::parameters::DRParameters* drParameters)
      : FrictionSolverDetails(drParameters) {}
  ~BaseFrictionSolver<Derived>() override = default;

  std::unique_ptr<FrictionSolver> clone() override {
    return std::make_unique<Derived>(*static_cast<Derived*>(this));
  }

  SEISSOL_DEVICE static void evaluatePoint(FrictionLawContext& ctx) {
    constexpr common::RangeType GpuRangeType{common::RangeType::GPU};

    const auto etaPDamp = ctx.data->drParameters.etaDampEnd > ctx.args->fullUpdateTime
                              ? ctx.data->drParameters.etaDamp
                              : 1.0;
    common::precomputeStressFromQInterpolated<GpuRangeType>(
        ctx.faultStresses,
        ctx.data->impAndEta[ctx.ltsFace],
        ctx.data->impedanceMatrices[ctx.ltsFace],
        ctx.data->qInterpolatedPlus[ctx.ltsFace],
        ctx.data->qInterpolatedMinus[ctx.ltsFace],
        etaPDamp,
        ctx.pointIndex);

    const auto isFrictionEnergyRequired{ctx.data->drParameters.isFrictionEnergyRequired};
    const auto isCheckAbortCriteraEnabled{ctx.data->drParameters.isCheckAbortCriteraEnabled};
    const auto devTerminatorSlipRateThreshold{ctx.data->drParameters.terminatorSlipRateThreshold};

    Derived::preHook(ctx);

    real startTime = 0;
    real updateTime = ctx.args->fullUpdateTime;
    for (uint32_t timeIndex = 0; timeIndex < misc::TimeSteps; ++timeIndex) {
      const real dt = ctx.args->deltaT[timeIndex];

      startTime = updateTime;
      updateTime += dt;

      for (uint32_t i = 0; i < ctx.data->drParameters.nucleationCount; ++i) {
        common::adjustInitialStress<GpuRangeType>(
            ctx.data->initialStressInFaultCS[ctx.ltsFace],
            ctx.data
                ->nucleationStressInFaultCS[ctx.ltsFace * ctx.data->drParameters.nucleationCount +
                                            i],
            ctx.data->initialPressure[ctx.ltsFace],
            ctx.data->nucleationPressure[ctx.ltsFace * ctx.data->drParameters.nucleationCount + i],
            updateTime,
            ctx.data->drParameters.t0[i],
            ctx.data->drParameters.s0[i],
            dt,
            ctx.pointIndex);
      }

      Derived::updateFrictionAndSlip(ctx, timeIndex);

      // time-dependent outputs
      common::saveRuptureFrontOutput<GpuRangeType>(ctx.data->ruptureTimePending[ctx.ltsFace],
                                                   ctx.data->ruptureTime[ctx.ltsFace],
                                                   ctx.data->slipRateMagnitude[ctx.ltsFace],
                                                   startTime,
                                                   ctx.pointIndex);

      Derived::saveDynamicStressOutput(ctx, startTime);

      common::savePeakSlipRateOutput<GpuRangeType>(ctx.data->slipRateMagnitude[ctx.ltsFace],
                                                   ctx.data->peakSlipRate[ctx.ltsFace],
                                                   ctx.pointIndex);

      if (isFrictionEnergyRequired && isCheckAbortCriteraEnabled) {
        common::updateTimeSinceSlipRateBelowThreshold<GpuRangeType>(
            ctx.data->slipRateMagnitude[ctx.ltsFace],
            ctx.data->ruptureTimePending[ctx.ltsFace],
            ctx.data->energyData[ctx.ltsFace],
            dt,
            devTerminatorSlipRateThreshold,
            ctx.pointIndex);
      }
    }

    Derived::postHook(ctx);

    common::postcomputeImposedStateFromNewStress<GpuRangeType>(
        ctx.faultStresses,
        ctx.tractionResults,
        ctx.data->impAndEta[ctx.ltsFace],
        ctx.data->impedanceMatrices[ctx.ltsFace],
        ctx.data->imposedStatePlus[ctx.ltsFace],
        ctx.data->imposedStateMinus[ctx.ltsFace],
        ctx.data->qInterpolatedPlus[ctx.ltsFace],
        ctx.data->qInterpolatedMinus[ctx.ltsFace],
        ctx.args->timeWeights,
        ctx.pointIndex);

    if (isFrictionEnergyRequired) {
      const auto energiesFromAcrossFaultVelocities{
          ctx.data->drParameters.energiesFromAcrossFaultVelocities};

      common::computeFrictionEnergy<GpuRangeType>(ctx.data->energyData[ctx.ltsFace],
                                                  ctx.data->qInterpolatedPlus[ctx.ltsFace],
                                                  ctx.data->qInterpolatedMinus[ctx.ltsFace],
                                                  ctx.data->impAndEta[ctx.ltsFace],
                                                  ctx.args->timeWeights,
                                                  ctx.args->spaceWeights,
                                                  ctx.data->godunovData[ctx.ltsFace],
                                                  ctx.data->slipRateMagnitude[ctx.ltsFace],
                                                  energiesFromAcrossFaultVelocities,
                                                  ctx.pointIndex);
    }
  }

  void setupLayer(DynamicRupture::Layer& layerData,
                  seissol::parallel::runtime::StreamRuntime& runtime) override {
    this->currLayerSize = layerData.size();
    FrictionSolverInterface::copyStorageToLocal(&dataHost, layerData);
    Derived::copySpecificStorageDataToLocal(&dataHost, layerData);
    dataHost.drParameters = *this->drParameters;
    device::DeviceInstance::getInstance().api->copyToAsync(
        data, &dataHost, sizeof(FrictionLawData), runtime.stream());
  }

  void evaluateKernel(seissol::parallel::runtime::StreamRuntime& runtime,
                      real fullUpdateTime,
                      const double* timeWeights,
                      const FrictionTime& frictionTime);

  void evaluate(real fullUpdateTime,
                const FrictionTime& frictionTime,
                const double* timeWeights,
                seissol::parallel::runtime::StreamRuntime& runtime) override {
    if (this->currLayerSize == 0) {
      return;
    }

    evaluateKernel(runtime, fullUpdateTime, timeWeights, frictionTime);
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_
