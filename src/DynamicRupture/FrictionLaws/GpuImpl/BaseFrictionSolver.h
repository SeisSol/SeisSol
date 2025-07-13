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
  int ltsFace;
  int pointIndex;
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
    constexpr common::RangeType gpuRangeType{common::RangeType::GPU};

    const auto etaPDamp =
        ctx.data->drParameters.etaStop > ctx.fullUpdateTime ? ctx.data->drParameters.etaHack : 1.0;
    common::precomputeStressFromQInterpolated<gpuRangeType>(
        ctx.faultStresses,
        ctx.data->impAndEta[ctx.ltsFace],
        ctx.data->impedanceMatrices[ctx.ltsFace],
        ctx.data->qInterpolatedPlus[ctx.ltsFace],
        ctx.data->qInterpolatedMinus[ctx.ltsFace],
        etaPDamp,
        ctx.pointIndex);

    Derived::preHook(ctx);

    real updateTime = ctx.fullUpdateTime;
    for (unsigned timeIndex = 0; timeIndex < ConvergenceOrder; ++timeIndex) {
      const real dt = ctx.data->deltaT[timeIndex];

      updateTime += dt;

      for (unsigned i = 0; i < ctx.data->drParameters.nucleationCount; ++i) {
        common::adjustInitialStress<gpuRangeType>(
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
    }
    Derived::postHook(ctx);

    common::saveRuptureFrontOutput<gpuRangeType>(ctx.data->ruptureTimePending[ctx.ltsFace],
                                                 ctx.data->ruptureTime[ctx.ltsFace],
                                                 ctx.data->slipRateMagnitude[ctx.ltsFace],
                                                 ctx.fullUpdateTime,
                                                 ctx.pointIndex);

    Derived::saveDynamicStressOutput(ctx);

    const auto devSumDt{ctx.data->sumDt};

    const auto isFrictionEnergyRequired{ctx.data->drParameters.isFrictionEnergyRequired};
    const auto isCheckAbortCriteraEnabled{ctx.data->drParameters.isCheckAbortCriteraEnabled};
    const auto devTerminatorSlipRateThreshold{ctx.data->drParameters.terminatorSlipRateThreshold};
    const auto energiesFromAcrossFaultVelocities{
        ctx.data->drParameters.energiesFromAcrossFaultVelocities};

    common::savePeakSlipRateOutput<gpuRangeType>(ctx.data->slipRateMagnitude[ctx.ltsFace],
                                                 ctx.data->peakSlipRate[ctx.ltsFace],
                                                 ctx.pointIndex);

    common::postcomputeImposedStateFromNewStress<gpuRangeType>(
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

      if (isCheckAbortCriteraEnabled) {
        common::updateTimeSinceSlipRateBelowThreshold<gpuRangeType>(
            ctx.data->slipRateMagnitude[ctx.ltsFace],
            ctx.data->ruptureTimePending[ctx.ltsFace],
            ctx.data->energyData[ctx.ltsFace],
            devSumDt,
            devTerminatorSlipRateThreshold,
            ctx.pointIndex);
      }

      common::computeFrictionEnergy<gpuRangeType>(ctx.data->energyData[ctx.ltsFace],
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
                      const double timeWeights[ConvergenceOrder]) {
    device::DeviceInstance::getInstance().api->copyToAsync(
        data, &dataHost, sizeof(FrictionLawData), runtime.stream());
    device::DeviceInstance::getInstance().api->copyToAsync(
        devTimeWeights, timeWeights, sizeof(double[ConvergenceOrder]), runtime.stream());
  }

  void evaluateKernel(seissol::parallel::runtime::StreamRuntime& runtime, real fullUpdateTime);

  void evaluate(seissol::initializer::Layer& layerData,
                const seissol::initializer::DynamicRupture* const dynRup,
                real fullUpdateTime,
                const double timeWeights[ConvergenceOrder],
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

    std::memcpy(dataHost.deltaT, deltaT, sizeof(decltype(deltaT)));
    dataHost.sumDt = sumDt;

    copyParameters(runtime, timeWeights);
    evaluateKernel(runtime, fullUpdateTime);
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_
