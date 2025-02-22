// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
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

namespace seissol::dr::friction_law::gpu {
struct InitialVariables {
  real absoluteShearTraction;
  real localSlipRate;
  real normalStress;
  real stateVarReference;
};

struct FrictionLawContext {
  int ltsFace;
  int pointIndex;
  FrictionLawData* __restrict data;

  FaultStresses<Executor::Device> faultStresses{};
  TractionResults<Executor::Device> tractionResults{};
  real fullUpdateTime;
  real stateVariableBuffer;
  real strengthBuffer;
  const double* __restrict devTimeWeights{nullptr};
  const real* __restrict devSpaceWeights{nullptr};
  const real* __restrict resampleMatrix{nullptr};
  real* __restrict deltaStateVar;
  InitialVariables initialVariables;

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

template <typename Derived> // , typename StdMath
class BaseFrictionSolver : public FrictionSolverDetails {
  public:
  explicit BaseFrictionSolver<Derived>(seissol::initializer::parameters::DRParameters* drParameters)
      : FrictionSolverDetails(drParameters) {}
  ~BaseFrictionSolver<Derived>() override = default;

  SEISSOL_DEVICE static void evaluatePoint(FrictionLawContext& ctx) {
    using StdMath = functions::HostStdFunctions;
    constexpr common::RangeType gpuRangeType{common::RangeType::GPU};

    auto* devImpAndEta{ctx.data->impAndEta};
    auto* devImpedanceMatrices{ctx.data->impedanceMatrices};
    auto* devQInterpolatedPlus{ctx.data->qInterpolatedPlus};
    auto* devQInterpolatedMinus{ctx.data->qInterpolatedMinus};

    common::precomputeStressFromQInterpolated<gpuRangeType>(ctx.faultStresses,
                                                            devImpAndEta[ctx.ltsFace],
                                                            devImpedanceMatrices[ctx.ltsFace],
                                                            devQInterpolatedPlus[ctx.ltsFace],
                                                            devQInterpolatedMinus[ctx.ltsFace],
                                                            ctx.pointIndex);

    Derived::preHook(ctx);
    for (unsigned timeIndex = 0; timeIndex < ConvergenceOrder; ++timeIndex) {
      const real t0{ctx.data->drParameters.t0};
      const real dt = ctx.data->deltaT[timeIndex];
      auto* devInitialStressInFaultCS{ctx.data->initialStressInFaultCS};
      const auto* devNucleationStressInFaultCS{ctx.data->nucleationStressInFaultCS};
      auto* devInitialPressure{ctx.data->initialPressure};
      const auto* devNucleationPressure{ctx.data->nucleationPressure};

      common::adjustInitialStress<gpuRangeType, StdMath>(devInitialStressInFaultCS[ctx.ltsFace],
                                                         devNucleationStressInFaultCS[ctx.ltsFace],
                                                         devInitialPressure[ctx.ltsFace],
                                                         devNucleationPressure[ctx.ltsFace],
                                                         ctx.fullUpdateTime,
                                                         t0,
                                                         dt,
                                                         ctx.pointIndex);

      Derived::updateFrictionAndSlip(ctx, timeIndex);
    }
    Derived::postHook(ctx);

    auto* devRuptureTimePending{ctx.data->ruptureTimePending};
    auto* devSlipRateMagnitude{ctx.data->slipRateMagnitude};
    auto* devRuptureTime{ctx.data->ruptureTime};

    common::saveRuptureFrontOutput<gpuRangeType>(devRuptureTimePending[ctx.ltsFace],
                                                 devRuptureTime[ctx.ltsFace],
                                                 devSlipRateMagnitude[ctx.ltsFace],
                                                 ctx.fullUpdateTime,
                                                 ctx.pointIndex);

    Derived::saveDynamicStressOutput(ctx);

    auto* devPeakSlipRate{ctx.data->peakSlipRate};
    auto* devImposedStatePlus{ctx.data->imposedStatePlus};
    auto* devImposedStateMinus{ctx.data->imposedStateMinus};
    auto* devEnergyData{ctx.data->energyData};
    auto* devGodunovData{ctx.data->godunovData};
    auto devSumDt{ctx.data->sumDt};

    auto isFrictionEnergyRequired{ctx.data->drParameters.isFrictionEnergyRequired};
    auto isCheckAbortCriteraEnabled{ctx.data->drParameters.isCheckAbortCriteraEnabled};
    auto devTerminatorSlipRateThreshold{ctx.data->drParameters.terminatorSlipRateThreshold};
    auto energiesFromAcrossFaultVelocities{
        ctx.data->drParameters.energiesFromAcrossFaultVelocities};

    common::savePeakSlipRateOutput<gpuRangeType>(
        devSlipRateMagnitude[ctx.ltsFace], devPeakSlipRate[ctx.ltsFace], ctx.pointIndex);

    common::postcomputeImposedStateFromNewStress<gpuRangeType>(ctx.faultStresses,
                                                               ctx.tractionResults,
                                                               devImpAndEta[ctx.ltsFace],
                                                               devImpedanceMatrices[ctx.ltsFace],
                                                               devImposedStatePlus[ctx.ltsFace],
                                                               devImposedStateMinus[ctx.ltsFace],
                                                               devQInterpolatedPlus[ctx.ltsFace],
                                                               devQInterpolatedMinus[ctx.ltsFace],
                                                               ctx.devTimeWeights,
                                                               ctx.pointIndex);

    if (isFrictionEnergyRequired) {

      if (isCheckAbortCriteraEnabled) {
        common::updateTimeSinceSlipRateBelowThreshold<gpuRangeType>(
            devSlipRateMagnitude[ctx.ltsFace],
            devRuptureTimePending[ctx.ltsFace],
            devEnergyData[ctx.ltsFace],
            devSumDt,
            devTerminatorSlipRateThreshold,
            ctx.pointIndex);
      }

      common::computeFrictionEnergy<gpuRangeType>(devEnergyData[ctx.ltsFace],
                                                  devQInterpolatedPlus[ctx.ltsFace],
                                                  devQInterpolatedMinus[ctx.ltsFace],
                                                  devImpAndEta[ctx.ltsFace],
                                                  ctx.devTimeWeights,
                                                  ctx.devSpaceWeights,
                                                  devGodunovData[ctx.ltsFace],
                                                  devSlipRateMagnitude[ctx.ltsFace],
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

    if (layerData.getNumberOfCells() == 0) {
      return;
    }

    // TODO: avoid copying the data all the time
    // TODO: allocate FrictionLawData as constant data

    FrictionSolverInterface::copyLtsTreeToLocal(&dataHost, layerData, dynRup, fullUpdateTime);
    Derived::copySpecificLtsDataTreeToLocal(&dataHost, layerData, dynRup, fullUpdateTime);
    this->currLayerSize = layerData.getNumberOfCells();
    dataHost.drParameters = *this->drParameters;

    std::memcpy(dataHost.deltaT, deltaT, sizeof(decltype(deltaT)));
    dataHost.sumDt = sumDt;

    copyParameters(runtime, timeWeights);
    evaluateKernel(runtime, fullUpdateTime);
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_
