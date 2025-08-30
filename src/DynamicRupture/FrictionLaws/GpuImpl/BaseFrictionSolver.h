// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_

#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverDetails.h"
#include "FrictionSolverInterface.h"
#include "Numerical/Functions.h"
#include <Common/Constants.h>
#include <DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h>
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

struct FrictionLawContext {
  std::size_t ltsFace;
  std::uint32_t pointIndex;
  FrictionLawData* __restrict data;

  FaultStresses<Executor::Device> faultStresses{};
  TractionResults<Executor::Device> tractionResults{};
  real fullUpdateTime{};
  real stateVariableBuffer{};
  real strengthBuffer{};
  const double* __restrict devTimeWeights{nullptr};
  const real* __restrict devSpaceWeights{nullptr};
  const real* __restrict resampleMatrix{nullptr};
  real* __restrict sharedMemory;
  InitialVariables initialVariables{};

  const real* __restrict TpInverseFourierCoefficients{nullptr};
  const real* __restrict TpGridPoints{nullptr};
  const real* __restrict HeatSource{nullptr};

  void* item;
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

  SEISSOL_DEVICE static void evaluatePoint(FrictionLawContext& ctx) {
    constexpr common::RangeType GpuRangeType{common::RangeType::GPU};

    const auto etaPDamp =
        ctx.data->drParameters.etaStop > ctx.fullUpdateTime ? ctx.data->drParameters.etaHack : 1.0;
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
    real updateTime = ctx.fullUpdateTime;
    for (uint32_t timeIndex = 0; timeIndex < misc::TimeSteps; ++timeIndex) {
      const real dt = ctx.data->deltaT[timeIndex];

      startTime = updateTime;
      updateTime += dt;

      for (uint32_t i = 0; i < ctx.data->drParameters.nucleationCount; ++i) {
        common::adjustInitialStress<GpuRangeType>(
            ctx.data->initialStressInFaultCS[ctx.ltsFace],
            ctx.data->nucleationStressInFaultCS[i][ctx.ltsFace],
            ctx.data->initialPressure[ctx.ltsFace],
            ctx.data->nucleationPressure[i][ctx.ltsFace],
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
        ctx.devTimeWeights,
        ctx.pointIndex);

    if (isFrictionEnergyRequired) {
      const auto energiesFromAcrossFaultVelocities{
          ctx.data->drParameters.energiesFromAcrossFaultVelocities};

      common::computeFrictionEnergy<GpuRangeType>(ctx.data->energyData[ctx.ltsFace],
                                                  ctx.data->qInterpolatedPlus[ctx.ltsFace],
                                                  ctx.data->qInterpolatedMinus[ctx.ltsFace],
                                                  ctx.data->impAndEta[ctx.ltsFace],
                                                  ctx.devTimeWeights,
                                                  ctx.devSpaceWeights,
                                                  ctx.data->godunovData[ctx.ltsFace],
                                                  ctx.data->slipRateMagnitude[ctx.ltsFace],
                                                  energiesFromAcrossFaultVelocities,
                                                  ctx.pointIndex);
    }
  }

  void copyParameters(seissol::parallel::runtime::StreamRuntime& runtime,
                      const double* timeWeights) {
    device::DeviceInstance::getInstance().api->copyToAsync(
        data, &dataHost, sizeof(FrictionLawData), runtime.stream());
    device::DeviceInstance::getInstance().api->copyToAsync(
        devTimeWeights, timeWeights, sizeof(double[misc::TimeSteps]), runtime.stream());
  }

  void evaluateKernel(seissol::parallel::runtime::StreamRuntime& runtime, real fullUpdateTime);

  void evaluate(seissol::initializer::Layer& layerData,
                const seissol::initializer::DynamicRupture* const dynRup,
                real fullUpdateTime,
                const FrictionTime& frictionTime,
                const double* timeWeights,
                seissol::parallel::runtime::StreamRuntime& runtime) override {

    if (layerData.size() == 0) {
      return;
    }

    // TODO: avoid copying the data all the time
    // TODO: allocate FrictionLawData as constant data

    FrictionSolverInterface::copyLtsTreeToLocal(&dataHost, layerData, dynRup, fullUpdateTime);
    Derived::copySpecificLtsDataTreeToLocal(&dataHost, layerData, dynRup, fullUpdateTime);
    this->currLayerSize = layerData.size();
    dataHost.drParameters = *this->drParameters;

    std::copy_n(frictionTime.deltaT.begin(), frictionTime.deltaT.size(), dataHost.deltaT);
    dataHost.sumDt = frictionTime.sumDt;

    copyParameters(runtime, timeWeights);
    evaluateKernel(runtime, fullUpdateTime);
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_
