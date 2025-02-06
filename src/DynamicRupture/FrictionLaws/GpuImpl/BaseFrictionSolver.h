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
#include "Numerical/SyclFunctions.h"
#include <Common/Constants.h>
#include <DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h>
#include <algorithm>

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
  FrictionLawData* data;

  FaultStresses<Executor::Device> faultStresses{};
  TractionResults<Executor::Device> tractionResults{};
  real stateVariableBuffer;
  real strengthBuffer;
  double* devTimeWeights{nullptr};
  real* devSpaceWeights{nullptr};
  real* resampleMatrix{nullptr};
  real* deltaStateVar;
  InitialVariables initialVariables;
  sycl::nd_item<1>* item;
};

template <typename Derived>
class BaseFrictionSolver : public FrictionSolverDetails {
  public:
  explicit BaseFrictionSolver<Derived>(seissol::initializer::parameters::DRParameters* drParameters)
      : FrictionSolverDetails(drParameters) {}
  ~BaseFrictionSolver<Derived>() override = default;

  void evaluate(seissol::initializer::Layer& layerData,
                const seissol::initializer::DynamicRupture* const dynRup,
                real fullUpdateTime,
                const double timeWeights[ConvergenceOrder],
                seissol::parallel::runtime::StreamRuntime& runtime) override {

    runtime.syncToSycl(&this->queue);

    // TODO: avoid copying the data all the time
    // TODO: allocate FrictionLawData as constant data

    FrictionSolverInterface::copyLtsTreeToLocal(&dataHost, layerData, dynRup, fullUpdateTime);
    Derived::copySpecificLtsDataTreeToLocal(&dataHost, layerData, dynRup, fullUpdateTime);
    this->currLayerSize = layerData.getNumberOfCells();
    dataHost.drParameters = *this->drParameters;

    std::memcpy(dataHost.deltaT, deltaT, sizeof(decltype(deltaT)));
    dataHost.sumDt = sumDt;

    this->queue.memcpy(data, &dataHost, sizeof(FrictionLawData));

    this->queue.memcpy(devTimeWeights, timeWeights, sizeof(double[ConvergenceOrder]));

    {
      constexpr common::RangeType gpuRangeType{common::RangeType::GPU};

      auto* data{this->data};
      auto* devTimeWeights{this->devTimeWeights};
      auto* devSpaceWeights{this->devSpaceWeights};
      auto* resampleMatrix{this->resampleMatrix};

      sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints}, {misc::NumPaddedPoints}};
      this->queue.submit([&](sycl::handler& cgh) {
        // NOLINTNEXTLINE
        sycl::accessor<real, 1, sycl::access::mode::read_write, sycl::access::target::local>
            deltaStateVar(misc::NumPaddedPoints, cgh);

        cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
          FrictionLawContext ctx{};
          ctx.data = data;
          ctx.devTimeWeights = devTimeWeights;
          ctx.devSpaceWeights = devSpaceWeights;
          ctx.resampleMatrix = resampleMatrix;
          ctx.deltaStateVar = &deltaStateVar[0];
          ctx.item = &item;

          auto* devImpAndEta{data->impAndEta};
          auto* devImpedanceMatrices{data->impedanceMatrices};
          auto* devQInterpolatedPlus{data->qInterpolatedPlus};
          auto* devQInterpolatedMinus{data->qInterpolatedMinus};

          const auto ltsFace = item.get_group().get_group_id(0);
          const auto pointIndex = item.get_local_id(0);

          ctx.ltsFace = ltsFace;
          ctx.pointIndex = pointIndex;

          common::precomputeStressFromQInterpolated<gpuRangeType>(ctx.faultStresses,
                                                                  devImpAndEta[ltsFace],
                                                                  devImpedanceMatrices[ltsFace],
                                                                  devQInterpolatedPlus[ltsFace],
                                                                  devQInterpolatedMinus[ltsFace],
                                                                  pointIndex);

          Derived::preHook(ctx);
          for (unsigned timeIndex = 0; timeIndex < ConvergenceOrder; ++timeIndex) {
            const real t0{data->drParameters.t0};
            const real dt = data->deltaT[timeIndex];
            auto* devInitialStressInFaultCS{data->initialStressInFaultCS};
            const auto* devNucleationStressInFaultCS{data->nucleationStressInFaultCS};
            auto* devInitialPressure{data->initialPressure};
            const auto* devNucleationPressure{data->nucleationPressure};

            using StdMath = seissol::functions::SyclStdFunctions;
            common::adjustInitialStress<gpuRangeType, StdMath>(
                devInitialStressInFaultCS[ltsFace],
                devNucleationStressInFaultCS[ltsFace],
                devInitialPressure[ltsFace],
                devNucleationPressure[ltsFace],
                fullUpdateTime,
                t0,
                dt,
                pointIndex);

            Derived::updateFrictionAndSlip(ctx, timeIndex);
          }
          Derived::postHook(ctx);

          auto* devRuptureTimePending{data->ruptureTimePending};
          auto* devSlipRateMagnitude{data->slipRateMagnitude};
          auto* devRuptureTime{data->ruptureTime};

          common::saveRuptureFrontOutput<gpuRangeType>(devRuptureTimePending[ltsFace],
                                                       devRuptureTime[ltsFace],
                                                       devSlipRateMagnitude[ltsFace],
                                                       fullUpdateTime,
                                                       pointIndex);

          Derived::saveDynamicStressOutput(ctx);

          auto* devPeakSlipRate{data->peakSlipRate};
          auto* devImposedStatePlus{data->imposedStatePlus};
          auto* devImposedStateMinus{data->imposedStateMinus};
          auto* devEnergyData{data->energyData};
          auto* devGodunovData{data->godunovData};
          auto devSumDt{data->sumDt};

          auto isFrictionEnergyRequired{data->drParameters.isFrictionEnergyRequired};
          auto isCheckAbortCriteraEnabled{data->drParameters.isCheckAbortCriteraEnabled};
          auto devTerminatorSlipRateThreshold{data->drParameters.terminatorSlipRateThreshold};

          common::savePeakSlipRateOutput<gpuRangeType>(
              devSlipRateMagnitude[ltsFace], devPeakSlipRate[ltsFace], pointIndex);

          common::postcomputeImposedStateFromNewStress<gpuRangeType>(ctx.faultStresses,
                                                                     ctx.tractionResults,
                                                                     devImpAndEta[ltsFace],
                                                                     devImpedanceMatrices[ltsFace],
                                                                     devImposedStatePlus[ltsFace],
                                                                     devImposedStateMinus[ltsFace],
                                                                     devQInterpolatedPlus[ltsFace],
                                                                     devQInterpolatedMinus[ltsFace],
                                                                     devTimeWeights,
                                                                     pointIndex);

          if (isFrictionEnergyRequired) {

            if (isCheckAbortCriteraEnabled) {
              common::updateTimeSinceSlipRateBelowThreshold<gpuRangeType>(
                  devSlipRateMagnitude[ltsFace],
                  devRuptureTimePending[ltsFace],
                  devEnergyData[ltsFace],
                  devSumDt,
                  devTerminatorSlipRateThreshold,
                  pointIndex);
            }

            common::computeFrictionEnergy<gpuRangeType>(devEnergyData[ltsFace],
                                                        devQInterpolatedPlus[ltsFace],
                                                        devQInterpolatedMinus[ltsFace],
                                                        devImpAndEta[ltsFace],
                                                        devTimeWeights,
                                                        devSpaceWeights,
                                                        devGodunovData[ltsFace],
                                                        pointIndex);
          }
        });
      });
    }

    runtime.syncFromSycl(&this->queue);
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_
