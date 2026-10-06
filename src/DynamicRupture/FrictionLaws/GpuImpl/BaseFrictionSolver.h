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
#include "Common/Real.h"
#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverDetails.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "FrictionSolverInterface.h"
#include "GeneratedCode/init.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Numerical/Functions.h"

#include <algorithm>

#ifdef SEISSOL_KERNELS_SYCL
#include <sycl/sycl.hpp>
#endif

namespace seissol::dr::friction_law::gpu {
template <typename Cfg>
struct InitialVariables {
  using real = Real<Cfg>;

  real absoluteShearTraction{};
  real localSlipRate{};
  real normalStress{};
  /// the same, before the slip rate dependent part and the clamp; the Newton solve needs it to
  /// evaluate sigma(V) itself
  real normalStressStick{};
  real stateVarReference{};
  real etaNormal{};
  /// unit slip direction; equals the normalized trial traction unless the impedance is anisotropic
  real slipDirection1{};
  real slipDirection2{};
};

template <typename Cfg>
struct FrictionLawArgs {
  using real = Real<Cfg>;

  const FrictionLawData<Cfg>* __restrict data{nullptr};
  const real* __restrict spaceWeights{nullptr};
  const real* __restrict resampleMatrix{nullptr};
  const real* __restrict tpInverseFourierCoefficients{nullptr};
  const real* __restrict tpGridPoints{nullptr};
  const real* __restrict heatSource{nullptr};

  real timeWeights[misc::TimeSteps<Cfg>]{};
  real deltaT[misc::TimeSteps<Cfg>]{};
  real fullUpdateTime{};
};

template <typename Cfg>
struct FrictionLawContext {
  using real = Real<Cfg>;

  std::size_t ltsFace{0};
  std::uint32_t pointIndex{0};
  const FrictionLawData<Cfg>* __restrict data{nullptr};
  const FrictionLawArgs<Cfg>* __restrict args{nullptr};

  real* __restrict sharedMemory{nullptr};
  void* item{nullptr};

  FaultStresses<Cfg, Executor::Device> faultStresses{};
  FaultStresses<Cfg, Executor::Device> initialStress{};
  TractionResults<Cfg, Executor::Device> tractionResults{};
  real stateVariableBuffer{};
  real strengthBuffer{};
  /// d(strength)/d(-sigma_eff); only used for the anisotropic normal/shear coupling
  real strengthSlopeBuffer{};
  InitialVariables<Cfg> initialVariables{};
};

#ifdef __CUDACC__
template <typename Cfg>
SEISSOL_DEVICE inline void deviceBarrier(FrictionLawContext<Cfg>& __restrict ctx) {
  __syncthreads();
}
template <typename Cfg>
SEISSOL_DEVICE inline void deviceWarpBarrier(FrictionLawContext<Cfg>& __restrict ctx) {
  __syncwarp();
}
template <typename Cfg>
SEISSOL_DEVICE inline bool deviceWarpAll(FrictionLawContext<Cfg>& __restrict ctx, bool value) {
  return __all_sync(0xffffffffU, static_cast<int>(value)) != 0;
}
#elif defined(__HIP__)
template <typename Cfg>
SEISSOL_DEVICE inline void deviceBarrier(FrictionLawContext<Cfg>& __restrict ctx) {
  __syncthreads();
}
template <typename Cfg>
SEISSOL_DEVICE inline void deviceWarpBarrier(FrictionLawContext<Cfg>& __restrict ctx) {
  // __syncwarp has no effect on current AMD GPUs (early 2026)
  // (nor does the HIP in our current CI support it)
}
template <typename Cfg>
SEISSOL_DEVICE inline bool deviceWarpAll(FrictionLawContext<Cfg>& __restrict ctx, bool value) {
  return __all(static_cast<int>(value)) != 0;
}
#elif defined(SEISSOL_KERNELS_SYCL)
template <typename Cfg>
inline void deviceBarrier(FrictionLawContext<Cfg>& __restrict ctx) {
  reinterpret_cast<sycl::nd_item<1>*>(ctx.item)->barrier(sycl::access::fence_space::local_space);
}
template <typename Cfg>
inline void deviceWarpBarrier(FrictionLawContext<Cfg>& __restrict ctx) {
  auto subgroup = reinterpret_cast<sycl::nd_item<1>*>(ctx.item)->get_sub_group();
  sycl::group_barrier(subgroup);
}
template <typename Cfg>
inline bool deviceWarpAll(FrictionLawContext<Cfg>& __restrict ctx, bool value) {
  auto subgroup = reinterpret_cast<sycl::nd_item<1>*>(ctx.item)->get_sub_group();
  return sycl::all_of_group(subgroup, value);
}
#else
template <typename Cfg>
inline void deviceBarrier(FrictionLawContext<Cfg>& __restrict /*ctx*/) {}
template <typename Cfg>
inline void deviceWarpBarrier(FrictionLawContext<Cfg>& __restrict /*ctx*/) {}
template <typename Cfg>
inline bool deviceWarpAll(FrictionLawContext<Cfg>& __restrict /*ctx*/, bool /*value*/) {
  return true;
}
#endif

template <typename Cfg>
SEISSOL_DEVICE inline Real<Cfg> resampleVariable(FrictionLawContext<Cfg>& __restrict ctx,
                                                 Real<Cfg> toResample) {
  constexpr auto Dim0 = misc::dimSize<init::resample<Cfg>, 0>();
  constexpr auto Dim1 = misc::dimSize<init::resample<Cfg>, 1>();
  static_assert(Dim0 == misc::NumPaddedPointsSingleSim<Cfg>);
  static_assert(Dim0 >= Dim1);

  ctx.sharedMemory[ctx.pointIndex] = toResample;
  deviceBarrier(ctx);

  const auto simPointIndex = ctx.pointIndex / Cfg::NumSimulations;
  const auto simId = ctx.pointIndex % Cfg::NumSimulations;
  constexpr uint32_t SimPointStride =
      multisim::MultisimHelperWrapper<Cfg>::MultisimEnabled ? Dim1 : 1U;
  constexpr uint32_t DataPointStride =
      multisim::MultisimHelperWrapper<Cfg>::MultisimEnabled ? 1U : Dim0;

  Real<Cfg> result{0};
  for (uint32_t i = 0; i < Dim1; ++i) {
    result += ctx.args->resampleMatrix[simPointIndex * SimPointStride + i * DataPointStride] *
              ctx.sharedMemory[i * Cfg::NumSimulations + simId];
  }
  deviceBarrier(ctx);

  return result;
}

template <typename Cfg, typename Derived>
class BaseFrictionSolver : public FrictionSolverDetails<Cfg> {
  public:
  using real = Real<Cfg>;

  explicit BaseFrictionSolver(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : FrictionSolverDetails<Cfg>(drParameters) {}
  ~BaseFrictionSolver() override = default;

  SEISSOL_DEVICE static void evaluatePoint(FrictionLawContext<Cfg>& __restrict ctx) {
    if constexpr (model::MaterialOf<Cfg>::SupportsDR) {
      constexpr common::RangeType GpuRangeType{common::RangeType::GPU};

      const auto etaPDamp = ctx.data->drParameters.etaDampEnd > ctx.args->fullUpdateTime
                                ? ctx.data->drParameters.etaDamp
                                : static_cast<real>(1.0);

      const auto isFrictionEnergyRequired{ctx.data->drParameters.isFrictionEnergyRequired};
      const auto isCheckAbortCriteraEnabled{ctx.data->drParameters.isCheckAbortCriteraEnabled};
      const auto devTerminatorSlipRateThreshold{ctx.data->drParameters.terminatorSlipRateThreshold};

      ImposedState<Cfg, Executor::Device> imposedState{};

      Derived::preHook(ctx);

      real startTime = 0;
      real updateTime = ctx.args->fullUpdateTime;

      for (uint32_t timeIndex = 0; timeIndex < misc::TimeSteps<Cfg>; ++timeIndex) {
        const real dt = ctx.args->deltaT[timeIndex];

        startTime = updateTime;
        updateTime += dt;

        common::precomputeStressFromQInterpolated<Cfg, GpuRangeType>(
            ctx.faultStresses,
            ctx.data->impAndEta[ctx.ltsFace],
            ctx.data->impedanceMatrices[ctx.ltsFace],
            ctx.data->qInterpolatedPlus[ctx.ltsFace],
            ctx.data->qInterpolatedMinus[ctx.ltsFace],
            etaPDamp,
            timeIndex,
            ctx.pointIndex);

        common::initializeTractionResults<Cfg, GpuRangeType>(
            ctx.faultStresses, ctx.tractionResults, ctx.pointIndex);

        const auto sourceCount = ctx.data->drParameters.sourceCount;
        common::computeInitialStress<Cfg, GpuRangeType>(
            ctx.initialStress,
            &ctx.data->stressSourceInFaultCS[ctx.ltsFace * sourceCount],
            &ctx.data->stressSourcePressure[ctx.ltsFace * sourceCount],
            &ctx.data->stressSourceOnset[ctx.ltsFace * sourceCount],
            &ctx.data->stressSourceRiseTime[ctx.ltsFace * sourceCount],
            sourceCount,
            updateTime,
            ctx.pointIndex);

        Derived::updateFrictionAndSlip(ctx, timeIndex);

        // time-dependent outputs
        common::saveRuptureFrontOutput<Cfg, GpuRangeType>(ctx.data->ruptureTimePending[ctx.ltsFace],
                                                          ctx.data->ruptureTime[ctx.ltsFace],
                                                          ctx.data->slipRateMagnitude[ctx.ltsFace],
                                                          startTime,
                                                          ctx.pointIndex);

        Derived::saveDynamicStressOutput(ctx, startTime);

        common::savePeakSlipRateOutput<Cfg, GpuRangeType>(ctx.data->slipRateMagnitude[ctx.ltsFace],
                                                          ctx.data->peakSlipRate[ctx.ltsFace],
                                                          ctx.pointIndex);

        if (isFrictionEnergyRequired && isCheckAbortCriteraEnabled) {
          common::updateTimeSinceSlipRateBelowThreshold<Cfg, GpuRangeType>(
              ctx.data->slipRateMagnitude[ctx.ltsFace],
              ctx.data->ruptureTimePending[ctx.ltsFace],
              ctx.data->energyData[ctx.ltsFace],
              dt,
              devTerminatorSlipRateThreshold,
              ctx.pointIndex);
        }

        common::postcomputeImposedStateFromNewStress<Cfg, GpuRangeType>(
            imposedState,
            ctx.faultStresses,
            ctx.tractionResults,
            ctx.data->impAndEta[ctx.ltsFace],
            ctx.data->impedanceMatrices[ctx.ltsFace],
            ctx.data->qInterpolatedPlus[ctx.ltsFace],
            ctx.data->qInterpolatedMinus[ctx.ltsFace],
            timeIndex,
            ctx.args->timeWeights[timeIndex],
            ctx.pointIndex);
      }

      Derived::postHook(ctx);

      common::finalizeImposedState<Cfg, GpuRangeType>(imposedState,
                                                      ctx.data->imposedStatePlus[ctx.ltsFace],
                                                      ctx.data->imposedStateMinus[ctx.ltsFace],
                                                      ctx.pointIndex);

      if (isFrictionEnergyRequired) {
        const auto energiesFromAcrossFaultVelocities{
            ctx.data->drParameters.energiesFromAcrossFaultVelocities};

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
  }

  void setupLayer(DynamicRupture::Layer& layerData,
                  seissol::parallel::runtime::StreamRuntime& runtime) override {
    this->currLayerSize_ = layerData.size();
    FrictionSolverInterface<Cfg>::copyStorageToLocal(&this->dataHost_, layerData);
    Derived::copySpecificStorageDataToLocal(&this->dataHost_, layerData);
    this->dataHost_.drParameters = this->drParameters_;
    device::DeviceInstance::instance().api().copyToAsync(
        this->data_, &this->dataHost_, sizeof(FrictionLawData<Cfg>), runtime.stream());
  }

  void evaluateKernel(seissol::parallel::runtime::StreamRuntime& runtime,
                      double fullUpdateTime,
                      const double* timeWeights,
                      const FrictionSolver::FrictionTime& frictionTime);

  void evaluate(double fullUpdateTime,
                const FrictionSolver::FrictionTime& frictionTime,
                const double* timeWeights,
                seissol::parallel::runtime::StreamRuntime& runtime) override {
    if (this->currLayerSize_ == 0) {
      return;
    }

    if constexpr (model::MaterialOf<Cfg>::SupportsDR) {
      evaluateKernel(runtime, fullUpdateTime, timeWeights, frictionTime);
    } else {
      logError() << "The material" << model::MaterialOf<Cfg>::Text
                 << "does not support DR friction law computations.";
    }
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_
