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
#include <Memory/Descriptor/DynamicRupture.h>
#include <algorithm>

#include "Common/Marker.h"

#ifdef SEISSOL_KERNELS_SYCL
#include <sycl/sycl.hpp>
#endif

namespace seissol::dr::friction_law::gpu {

template<typename RealT>
struct InitialVariables {
  RealT absoluteShearTraction{};
  RealT localSlipRate{};
  RealT normalStress{};
  RealT stateVarReference{};
};

template<typename Cfg>
struct FrictionLawArgs {
  const FrictionLawData* __restrict data{nullptr};
  const Real<Cfg>* __restrict spaceWeights{nullptr};
  const Real<Cfg>* __restrict resampleMatrix{nullptr};
  const Real<Cfg>* __restrict tpInverseFourierCoefficients{nullptr};
  const Real<Cfg>* __restrict tpGridPoints{nullptr};
  const Real<Cfg>* __restrict heatSource{nullptr};

  Real<Cfg> fullUpdateTime;
  double timeWeights[Cfg::ConvergenceOrder];
  Real<Cfg> deltaT[Cfg::ConvergenceOrder];
  Real<Cfg> sumDt;
};

template<typename Cfg>
struct FrictionLawContext {
  std::size_t ltsFace;
  std::uint32_t pointIndex;
  const FrictionLawData* __restrict data;
  const FrictionLawArgs<Cfg>* __restrict args;

  Real<Cfg>* __restrict sharedMemory;
  void* item;

  FaultStresses<Cfg, Executor::Device> faultStresses{};
  TractionResults<Cfg, Executor::Device> tractionResults{};
  Real<Cfg> stateVariableBuffer{};
  Real<Cfg> strengthBuffer{};
  InitialVariables<Real<Cfg>> initialVariables{};
};

#ifdef __CUDACC__
template<typename Cfg>
SEISSOL_DEVICE inline void deviceBarrier(FrictionLawContext<Cfg>& ctx) { __syncthreads(); }
#elif defined(__HIP__)
template<typename Cfg>
SEISSOL_DEVICE inline void deviceBarrier(FrictionLawContext<Cfg>& ctx) { __syncthreads(); }
#elif defined(SEISSOL_KERNELS_SYCL)
template<typename Cfg>
inline void deviceBarrier(FrictionLawContext<Cfg>& ctx) {
  reinterpret_cast<sycl::nd_item<1>*>(ctx.item)->barrier(sycl::access::fence_space::local_space);
}
#else
template<typename Cfg>
inline void deviceBarrier(FrictionLawContext<Cfg>& ctx) {}
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

  SEISSOL_DEVICE static void evaluatePoint(FrictionLawContext<Cfg>& ctx) {
    constexpr common::RangeType GpuRangeType{common::RangeType::GPU};

    const auto etaPDamp = ctx.data->drParameters.etaStop > ctx.args->fullUpdateTime
                              ? ctx.data->drParameters.etaHack
                              : 1.0;
    common::precomputeStressFromQInterpolated<Cfg, GpuRangeType>(
        ctx.faultStresses,
        ctx.data->impAndEta[ctx.ltsFace],
        ctx.data->impedanceMatrices[ctx.ltsFace],
        ctx.data->qInterpolatedPlus[ctx.ltsFace],
        ctx.data->qInterpolatedMinus[ctx.ltsFace],
        etaPDamp,
        ctx.pointIndex);

    Derived::preHook(ctx);

    real updateTime = ctx.args->fullUpdateTime;
    for (uint32_t timeIndex = 0; timeIndex < Cfg::ConvergenceOrder; ++timeIndex) {
      const real dt = ctx.args->deltaT[timeIndex];

      updateTime += dt;

      for (uint32_t i = 0; i < ctx.data->drParameters.nucleationCount; ++i) {
        common::adjustInitialStress<Cfg, GpuRangeType>(
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
    }
    Derived::postHook(ctx);

    common::saveRuptureFrontOutput<Cfg, GpuRangeType>(ctx.data->ruptureTimePending[ctx.ltsFace],
                                                 ctx.data->ruptureTime[ctx.ltsFace],
                                                 ctx.data->slipRateMagnitude[ctx.ltsFace],
                                                 ctx.args->fullUpdateTime,
                                                 ctx.pointIndex);

    Derived::saveDynamicStressOutput(ctx);

    const auto devSumDt{ctx.args->sumDt};

    const auto isFrictionEnergyRequired{ctx.data->drParameters.isFrictionEnergyRequired};
    const auto isCheckAbortCriteraEnabled{ctx.data->drParameters.isCheckAbortCriteraEnabled};
    const auto devTerminatorSlipRateThreshold{ctx.data->drParameters.terminatorSlipRateThreshold};
    const auto energiesFromAcrossFaultVelocities{
        ctx.data->drParameters.energiesFromAcrossFaultVelocities};

    common::savePeakSlipRateOutput<Cfg, GpuRangeType>(ctx.data->slipRateMagnitude[ctx.ltsFace],
                                                 ctx.data->peakSlipRate[ctx.ltsFace],
                                                 ctx.pointIndex);

    common::postcomputeImposedStateFromNewStress<Cfg, GpuRangeType>(
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

      if (isCheckAbortCriteraEnabled) {
        common::updateTimeSinceSlipRateBelowThreshold<Cfg, GpuRangeType>(
            ctx.data->slipRateMagnitude[ctx.ltsFace],
            ctx.data->ruptureTimePending[ctx.ltsFace],
            ctx.data->energyData[ctx.ltsFace],
            devSumDt,
            devTerminatorSlipRateThreshold,
            ctx.pointIndex);
      }

      common::computeFrictionEnergy<Cfg, GpuRangeType>(ctx.data->energyData[ctx.ltsFace],
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
                      double fullUpdateTime,
                      const double timeWeights[Cfg::ConvergenceOrder],
                      const FrictionTime& frictionTime);

  void evaluate(double fullUpdateTime,
                const FrictionTime& frictionTime,
                const double timeWeights[Cfg::ConvergenceOrder],
                seissol::parallel::runtime::StreamRuntime& runtime) override {
    if (this->currLayerSize == 0) {
      return;
    }

    evaluateKernel(runtime, fullUpdateTime, timeWeights, frictionTime);
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_
