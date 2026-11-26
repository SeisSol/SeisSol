// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_FRICTIONLAWDRIVER_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_FRICTIONLAWDRIVER_H_

#include "Core/Common/Marker.h"
#include "Core/Common/Real.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawContext.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolverCommon.h"
#include "Features/DynamicRupture/Typedefs.h"
#include "WavePropagation/Physics/Datastructures.h"

#include <cstdint>

namespace seissol::dr::friction_law {

/**
 * A time step of the friction law `LawT` at the points of a range of a face, on the host and the
 * device alike (see FrictionLawContext): for each sub-step, the stresses from the interpolated
 * unknowns, the friction and the slip (LawT::updateFrictionAndSlip), the outputs and the imposed
 * state; then the energies. The law supplies the hooks:
 *   SlipRateAtIntervalEnds, preHook(ctx), updateFrictionAndSlip(ctx, timeIndex),
 *   dynamicStressQuantity(ctx), saveDynamicStressOutput(ctx, previous, times), postHook(ctx).
 */
template <typename Cfg, typename LawT>
struct FrictionLawDriver {
  using real = Real<Cfg>;

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void evaluate(FrictionLawContext<Cfg, Type>& __restrict ctx) {
    if constexpr (model::MaterialOf<Cfg>::SupportsDR) {
      const auto etaPDamp = ctx.data->drParameters.etaDampEnd > ctx.fullUpdateTime
                                ? ctx.data->drParameters.etaDamp
                                : static_cast<real>(1.0);

      const auto isFrictionEnergyRequired{ctx.data->drParameters.isFrictionEnergyRequired};
      const auto isCheckAbortCriteraEnabled{ctx.data->drParameters.isCheckAbortCriteraEnabled};
      const auto terminatorSlipRateThreshold{ctx.data->drParameters.terminatorSlipRateThreshold};

      alignas(pointAlignment<Cfg, Type>()) ImposedState<Cfg, FrictionLawContext<Cfg, Type>::Exec>
          imposedState{};

      LawT::preHook(ctx);

      real updateTime = ctx.fullUpdateTime;

      for (uint32_t timeIndex = 0; timeIndex < misc::TimeSteps<Cfg>; ++timeIndex) {
        const real dt = ctx.args->deltaT[timeIndex];

        updateTime += dt;
        const auto times =
            common::subStepTimes<Cfg>(ctx.fullUpdateTime, ctx.args->deltaT, timeIndex);
        const auto rateTimes = common::slipRateTimes<Cfg>(times, LawT::SlipRateAtIntervalEnds);
        // the values before the update, for the times of the events that it crosses the threshold
        // of (see common::crossingTime)
        alignas(pointAlignment<Cfg, Type>()) PointValues<Cfg, Type> previousSlipRate{};
        forEachPoint(ctx, [&](std::uint32_t pointIndex) {
          previousSlipRate[pointIndex] = ctx.data->slipRateMagnitude[ctx.ltsFace][pointIndex];
        });
        alignas(pointAlignment<Cfg, Type>()) const PointValues<Cfg, Type>
            previousDynamicStressQuantity = LawT::dynamicStressQuantity(ctx);

        common::precomputeStressFromQInterpolated<Cfg, Type>(
            ctx.faultStresses,
            ctx.data->impAndEta[ctx.ltsFace],
            ctx.data->impedanceMatrices[ctx.ltsFace],
            ctx.data->qInterpolatedPlus[ctx.ltsFace],
            ctx.data->qInterpolatedMinus[ctx.ltsFace],
            etaPDamp,
            timeIndex,
            firstPoint(ctx));

        common::initializeTractionResults<Cfg, Type>(
            ctx.faultStresses, ctx.tractionResults, firstPoint(ctx));

        const auto sourceCount = ctx.data->drParameters.sourceCount;
        common::computeInitialStress<Cfg, Type>(
            ctx.initialStress,
            &ctx.data->stressSourceInFaultCS[ctx.ltsFace * sourceCount],
            &ctx.data->stressSourcePressure[ctx.ltsFace * sourceCount],
            &ctx.data->stressSourceOnset[ctx.ltsFace * sourceCount],
            &ctx.data->stressSourceRiseTime[ctx.ltsFace * sourceCount],
            sourceCount,
            updateTime,
            firstPoint(ctx));

        LawT::updateFrictionAndSlip(ctx, timeIndex);

        // time-dependent outputs
        common::saveRuptureFrontOutput<Cfg, Type>(ctx.data->ruptureTimePending[ctx.ltsFace],
                                                  ctx.data->ruptureTime[ctx.ltsFace],
                                                  &previousSlipRate[firstPoint(ctx)],
                                                  ctx.data->slipRateMagnitude[ctx.ltsFace],
                                                  rateTimes.previousPoint,
                                                  rateTimes.point,
                                                  firstPoint(ctx));

        LawT::saveDynamicStressOutput(ctx, previousDynamicStressQuantity, times);

        common::savePeakSlipRateOutput<Cfg, Type>(ctx.data->slipRateMagnitude[ctx.ltsFace],
                                                  ctx.data->peakSlipRate[ctx.ltsFace],
                                                  firstPoint(ctx));

        if (isFrictionEnergyRequired && isCheckAbortCriteraEnabled) {
          common::updateTimeSinceSlipRateBelowThreshold<Cfg, Type>(
              ctx.data->slipRateMagnitude[ctx.ltsFace],
              ctx.data->ruptureTimePending[ctx.ltsFace],
              ctx.data->energyData[ctx.ltsFace],
              dt,
              terminatorSlipRateThreshold,
              firstPoint(ctx));
        }

        common::postcomputeImposedStateFromNewStress<Cfg, Type>(
            imposedState,
            ctx.faultStresses,
            ctx.tractionResults,
            ctx.data->impAndEta[ctx.ltsFace],
            ctx.data->impedanceMatrices[ctx.ltsFace],
            ctx.data->qInterpolatedPlus[ctx.ltsFace],
            ctx.data->qInterpolatedMinus[ctx.ltsFace],
            timeIndex,
            ctx.args->timeWeights[timeIndex],
            firstPoint(ctx));
      }

      LawT::postHook(ctx);

      common::finalizeImposedState<Cfg, Type>(imposedState,
                                              ctx.data->imposedStatePlus[ctx.ltsFace],
                                              ctx.data->imposedStateMinus[ctx.ltsFace],
                                              firstPoint(ctx));

      if (isFrictionEnergyRequired) {
        const auto energiesFromAcrossFaultVelocities{
            ctx.data->drParameters.energiesFromAcrossFaultVelocities};

        common::computeFrictionEnergy<Cfg, Type>(ctx.data->energyData[ctx.ltsFace],
                                                 ctx.data->qInterpolatedPlus[ctx.ltsFace],
                                                 ctx.data->qInterpolatedMinus[ctx.ltsFace],
                                                 ctx.data->impAndEta[ctx.ltsFace],
                                                 ctx.args->timeWeights,
                                                 ctx.args->spaceWeights,
                                                 ctx.data->godunovData[ctx.ltsFace],
                                                 ctx.data->slipRateMagnitude[ctx.ltsFace],
                                                 energiesFromAcrossFaultVelocities,
                                                 firstPoint(ctx));
      }
    }
  }
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_FRICTIONLAWDRIVER_H_
