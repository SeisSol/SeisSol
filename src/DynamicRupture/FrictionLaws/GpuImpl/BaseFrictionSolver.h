// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_

#include "Common/Constants.h"
#include "Common/Marker.h"
#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverDetails.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "FrictionSolverInterface.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Numerical/Functions.h"

#include <algorithm>

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

  real fullUpdateTime{};
  double timeWeights[misc::TimeSteps]{};
  real deltaT[misc::TimeSteps]{};
};

struct FrictionLawContext {
  std::size_t ltsFace{0};
  std::uint32_t pointIndex{0};
  const FrictionLawData* __restrict data{nullptr};
  const FrictionLawArgs* __restrict args{nullptr};

  real* __restrict sharedMemory{nullptr};
  void* item{nullptr};

  FaultStresses<Executor::Device> faultStresses{};
  TractionResults<Executor::Device> tractionResults{};
  real stateVariableBuffer{};
  real strengthBuffer{};
  InitialVariables initialVariables{};
};

#ifdef __CUDACC__
SEISSOL_DEVICE inline void deviceBarrier(FrictionLawContext& /*ctx*/) { __syncthreads(); }
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
  explicit BaseFrictionSolver(seissol::initializer::parameters::DRParameters* drParameters)
      : FrictionSolverDetails(drParameters) {}
  ~BaseFrictionSolver() override = default;

  std::unique_ptr<FrictionSolver> clone() override {
    return std::make_unique<Derived>(*static_cast<Derived*>(this));
  }

  SEISSOL_DEVICE static void evaluatePoint(FrictionLawContext& ctx) {
    if constexpr (model::MaterialT::SupportsDR) {
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
              ctx.data
                  ->nucleationPressure[ctx.ltsFace * ctx.data->drParameters.nucleationCount + i],
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
  }

  void setupLayer(DynamicRupture::Layer& layerData,
                  seissol::parallel::runtime::StreamRuntime& runtime) override {
    this->currLayerSize_ = layerData.size();
    FrictionSolverInterface::copyStorageToLocal(&dataHost_, layerData);
    Derived::copySpecificStorageDataToLocal(&dataHost_, layerData);
    dataHost_.drParameters = *this->drParameters_;
    device::DeviceInstance::getInstance().api->copyToAsync(
        data_, &dataHost_, sizeof(FrictionLawData), runtime.stream());
  }

  void evaluateKernel(seissol::parallel::runtime::StreamRuntime& runtime,
                      real fullUpdateTime,
                      const double* timeWeights,
                      const FrictionTime& frictionTime);

  void evaluate(real fullUpdateTime,
                const FrictionTime& frictionTime,
                const double* timeWeights,
                seissol::parallel::runtime::StreamRuntime& runtime) override {
    if (this->currLayerSize_ == 0) {
      return;
    }

    if constexpr (model::MaterialT::SupportsDR) {
      evaluateKernel(runtime, fullUpdateTime, timeWeights, frictionTime);
    } else {
      logError() << "The material" << model::MaterialT::Text
                 << "does not support DR friction law computations.";
    }
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_
