// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_LINEARSLIPWEAKENING_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_LINEARSLIPWEAKENING_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include <DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h>

namespace seissol::dr::friction_law::gpu {

/**
 * Abstract Class implementing the general structure of linear slip weakening friction laws.
 * specific implementation is done by overriding and implementing the hook functions (via CRTP).
 */
template <typename Derived>
class LinearSlipWeakeningBase : public BaseFrictionSolver<LinearSlipWeakeningBase<Derived>> {
  public:
  LinearSlipWeakeningBase<Derived>(seissol::initializer::parameters::DRParameters* drParameters)
      : BaseFrictionSolver<LinearSlipWeakeningBase<Derived>>(drParameters){};

  void allocateAuxiliaryMemory() override { FrictionSolverDetails::allocateAuxiliaryMemory(); }

  static void
      copySpecificLtsDataTreeToLocal(FrictionLawData* data,
                                     seissol::initializer::Layer& layerData,
                                     const seissol::initializer::DynamicRupture* const dynRup,
                                     real fullUpdateTime) {
    Derived::copySpecificLtsDataTreeToLocal(data, layerData, dynRup, fullUpdateTime);
  }

  SEISSOL_DEVICE static void updateFrictionAndSlip(FrictionLawContext& ctx, unsigned timeIndex) {
    // computes fault strength, which is the critical value whether active slip exists.
    Derived::calcStrengthHook(ctx, timeIndex);
    // computes resulting slip rates, traction and slip dependent on current friction
    // coefficient and strength
    calcSlipRateAndTraction(ctx, timeIndex);
    Derived::calcStateVariableHook(ctx, timeIndex);
    frictionFunctionHook(ctx);
  }

  /**
   *  compute the slip rate and the traction from the fault strength and fault stresses
   *  also updates the directional slip1 and slip2
   */
  SEISSOL_DEVICE static void calcSlipRateAndTraction(FrictionLawContext& ctx,
                                                     unsigned int timeIndex) {
    const auto& devImpAndEta{ctx.data->impAndEta[ctx.ltsFace]};
    const auto deltaT{ctx.data->deltaT[timeIndex]};

    auto& faultStresses = ctx.faultStresses;
    auto& tractionResults = ctx.tractionResults;
    auto& strength = ctx.strengthBuffer;

    // calculate absolute value of stress in Y and Z direction
    const real totalStress1 = ctx.data->initialStressInFaultCS[ctx.ltsFace][ctx.pointIndex][3] +
                              faultStresses.traction1[timeIndex];
    const real totalStress2 = ctx.data->initialStressInFaultCS[ctx.ltsFace][ctx.pointIndex][5] +
                              faultStresses.traction2[timeIndex];
    const real absoluteShearStress = misc::magnitude(totalStress1, totalStress2);
    // calculate slip rates
    ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex] =
        std::max(static_cast<real>(0.0), (absoluteShearStress - strength) * devImpAndEta.invEtaS);
    const auto divisor =
        strength + devImpAndEta.etaS * ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex];
    ctx.data->slipRate1[ctx.ltsFace][ctx.pointIndex] =
        ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex] * totalStress1 / divisor;
    ctx.data->slipRate2[ctx.ltsFace][ctx.pointIndex] =
        ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex] * totalStress2 / divisor;
    // calculate traction
    tractionResults.traction1[timeIndex] =
        faultStresses.traction1[timeIndex] -
        devImpAndEta.etaS * ctx.data->slipRate1[ctx.ltsFace][ctx.pointIndex];
    tractionResults.traction2[timeIndex] =
        faultStresses.traction2[timeIndex] -
        devImpAndEta.etaS * ctx.data->slipRate2[ctx.ltsFace][ctx.pointIndex];
    ctx.data->traction1[ctx.ltsFace][ctx.pointIndex] = tractionResults.traction1[timeIndex];
    ctx.data->traction2[ctx.ltsFace][ctx.pointIndex] = tractionResults.traction2[timeIndex];
    // update directional slip
    ctx.data->slip1[ctx.ltsFace][ctx.pointIndex] +=
        ctx.data->slipRate1[ctx.ltsFace][ctx.pointIndex] * deltaT;
    ctx.data->slip2[ctx.ltsFace][ctx.pointIndex] +=
        ctx.data->slipRate2[ctx.ltsFace][ctx.pointIndex] * deltaT;
  }

  /**
   * evaluate friction law: updated mu -> friction law
   * for example see Carsten Uphoff's thesis: Eq. 2.45
   */
  SEISSOL_DEVICE static void frictionFunctionHook(FrictionLawContext& ctx) {
    auto& stateVariable = ctx.stateVariableBuffer;
    ctx.data->mu[ctx.ltsFace][ctx.pointIndex] =
        ctx.data->muS[ctx.ltsFace][ctx.pointIndex] -
        (ctx.data->muS[ctx.ltsFace][ctx.pointIndex] - ctx.data->muD[ctx.ltsFace][ctx.pointIndex]) *
            stateVariable;
    // instantaneous healing
    if ((ctx.data->peakSlipRate[ctx.ltsFace][ctx.pointIndex] >
         ctx.data->drParameters.healingThreshold) &&
        (ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex] <
         ctx.data->drParameters.healingThreshold)) {
      ctx.data->mu[ctx.ltsFace][ctx.pointIndex] = ctx.data->muS[ctx.ltsFace][ctx.pointIndex];
      stateVariable = 0.0;
    }
  }

  /*
   * output time when shear stress is equal to the dynamic stress after rupture arrived
   * currently only for linear slip weakening
   */
  SEISSOL_DEVICE static void saveDynamicStressOutput(FrictionLawContext& ctx) {
    const auto fullUpdateTime{ctx.data->mFullUpdateTime};

    if (ctx.data->dynStressTimePending[ctx.ltsFace][ctx.pointIndex] &&
        std::fabs(ctx.data->accumulatedSlipMagnitude[ctx.ltsFace][ctx.pointIndex]) >=
            ctx.data->dC[ctx.ltsFace][ctx.pointIndex]) {
      ctx.data->dynStressTime[ctx.ltsFace][ctx.pointIndex] = fullUpdateTime;
      ctx.data->dynStressTimePending[ctx.ltsFace][ctx.pointIndex] = false;
    }
  }

  SEISSOL_DEVICE static void preHook(FrictionLawContext& ctx) {}
  SEISSOL_DEVICE static void postHook(FrictionLawContext& ctx) {}

  protected:
  static constexpr real U0 = 10e-14;
};

template <class SpecializationT>
class LinearSlipWeakeningLaw
    : public LinearSlipWeakeningBase<LinearSlipWeakeningLaw<SpecializationT>> {
  public:
  LinearSlipWeakeningLaw<SpecializationT>(
      seissol::initializer::parameters::DRParameters* drParameters)
      : LinearSlipWeakeningBase<LinearSlipWeakeningLaw<SpecializationT>>(drParameters),
        specialization(drParameters){};

  static void
      copySpecificLtsDataTreeToLocal(FrictionLawData* data,
                                     seissol::initializer::Layer& layerData,
                                     const seissol::initializer::DynamicRupture* const dynRup,
                                     real fullUpdateTime) {
    const auto* concreteLts =
        dynamic_cast<const seissol::initializer::LTSLinearSlipWeakening*>(dynRup);
    data->dC = layerData.var(concreteLts->dC, seissol::initializer::AllocationPlace::Device);
    data->muS = layerData.var(concreteLts->muS, seissol::initializer::AllocationPlace::Device);
    data->muD = layerData.var(concreteLts->muD, seissol::initializer::AllocationPlace::Device);
    data->cohesion =
        layerData.var(concreteLts->cohesion, seissol::initializer::AllocationPlace::Device);
    data->forcedRuptureTime = layerData.var(concreteLts->forcedRuptureTime,
                                            seissol::initializer::AllocationPlace::Device);
    SpecializationT::copyLtsTreeToLocal(data, layerData, dynRup, fullUpdateTime);
  }

  SEISSOL_DEVICE static void calcStrengthHook(FrictionLawContext& ctx, unsigned int timeIndex) {

    const auto deltaT{ctx.data->deltaT[timeIndex]};

    const auto vStar{ctx.data->drParameters.vStar};
    const auto prakashLength{ctx.data->drParameters.prakashLength};

    auto& faultStresses = ctx.faultStresses;
    auto& strength = ctx.strengthBuffer;

    const real totalNormalStress =
        ctx.data->initialStressInFaultCS[ctx.ltsFace][ctx.pointIndex][0] +
        faultStresses.normalStress[timeIndex];
    strength = -ctx.data->cohesion[ctx.ltsFace][ctx.pointIndex] -
               ctx.data->mu[ctx.ltsFace][ctx.pointIndex] *
                   std::min(totalNormalStress, static_cast<real>(0.0));

    strength =
        SpecializationT::strengthHook(ctx,
                                      strength,
                                      ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex],
                                      deltaT,
                                      vStar,
                                      prakashLength);
  }

  SEISSOL_DEVICE static void calcStateVariableHook(FrictionLawContext& ctx,
                                                   unsigned int timeIndex) {
    const auto deltaT{ctx.data->deltaT[timeIndex]};
    const real tn{ctx.data->mFullUpdateTime + deltaT};
    const auto t0{ctx.data->drParameters.t0};
    const auto tpProxyExponent{ctx.data->drParameters.tpProxyExponent};

    const real resampledSlipRate =
        SpecializationT::resampleSlipRate(ctx, ctx.data->slipRateMagnitude[ctx.ltsFace]);

    // integrate slip rate to get slip = state variable
    ctx.data->accumulatedSlipMagnitude[ctx.ltsFace][ctx.pointIndex] += resampledSlipRate * deltaT;

    // Actually slip is already the stateVariable for this FL, but to simplify the next
    // equations we divide it here by the critical distance.
    const real localStateVariable = SpecializationT::stateVariableHook(
        ctx,
        ctx.data->accumulatedSlipMagnitude[ctx.ltsFace][ctx.pointIndex],
        ctx.data->dC[ctx.ltsFace][ctx.pointIndex],
        tpProxyExponent);

    real f2 = 0.0;
    if (t0 == 0) {
      f2 =
          1.0 * static_cast<double>(tn >= ctx.data->forcedRuptureTime[ctx.ltsFace][ctx.pointIndex]);
    } else {
      f2 = misc::clamp((tn - ctx.data->forcedRuptureTime[ctx.ltsFace][ctx.pointIndex]) / t0,
                       static_cast<real>(0.0),
                       static_cast<real>(1.0));
    }
    ctx.stateVariableBuffer = std::max(localStateVariable, f2);
  }

  protected:
  SpecializationT specialization;
};

class NoSpecialization {
  public:
  NoSpecialization(seissol::initializer::parameters::DRParameters* parameters) {};

  static void copyLtsTreeToLocal(FrictionLawData* data,
                                 seissol::initializer::Layer& layerData,
                                 const seissol::initializer::DynamicRupture* const dynRup,
                                 real fullUpdateTime) {}

  SEISSOL_DEVICE static real
      resampleSlipRate(FrictionLawContext& ctx,
                       const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints]) {

    // perform matrix vector multiplication

    constexpr auto Dim0 = misc::dimSize<init::resample, 0>();
    constexpr auto Dim1 = misc::dimSize<init::resample, 1>();
    static_assert(Dim0 == misc::NumPaddedPoints);
    static_assert(Dim0 >= Dim1);

    ctx.sharedMemory[ctx.pointIndex] = slipRateMagnitude[ctx.pointIndex];
    deviceBarrier(ctx);

    real result{0.0};
    for (size_t i{0}; i < Dim1; ++i) {
      result += ctx.resampleMatrix[ctx.pointIndex + i * Dim0] * ctx.sharedMemory[i];
    }
    return result;
  };

  SEISSOL_DEVICE static real stateVariableHook(FrictionLawContext& ctx,
                                               real localAccumulatedSlip,
                                               real localDc,
                                               real tpProxyExponent) {
    return std::min(std::fabs(localAccumulatedSlip) / localDc, static_cast<real>(1.0));
  };

  SEISSOL_DEVICE static real strengthHook(FrictionLawContext& ctx,
                                          real strength,
                                          real localSlipRate,
                                          real deltaT,
                                          real vStar,
                                          real prakashLength) {
    return strength;
  };
};

class BiMaterialFault {
  public:
  BiMaterialFault(seissol::initializer::parameters::DRParameters* parameters) {};

  static void copyLtsTreeToLocal(FrictionLawData* data,
                                 seissol::initializer::Layer& layerData,
                                 const seissol::initializer::DynamicRupture* const dynRup,
                                 real fullUpdateTime) {
    const auto* concreteLts =
        dynamic_cast<const seissol::initializer::LTSLinearSlipWeakeningBimaterial*>(dynRup);
    data->regularizedStrength = layerData.var(concreteLts->regularizedStrength,
                                              seissol::initializer::AllocationPlace::Device);
  }

  SEISSOL_DEVICE static real
      resampleSlipRate(FrictionLawContext& ctx,
                       const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints]) {
    return slipRateMagnitude[ctx.pointIndex];
  };

  SEISSOL_DEVICE static real stateVariableHook(FrictionLawContext& ctx,
                                               real localAccumulatedSlip,
                                               real localDc,
                                               real tpProxyExponent) {
    return std::min(std::fabs(localAccumulatedSlip) / localDc, static_cast<real>(1.0));
  };

  SEISSOL_DEVICE static real strengthHook(FrictionLawContext& ctx,
                                          real faultStrength,
                                          real localSlipRate,
                                          real deltaT,
                                          real vStar,
                                          real prakashLength) {
    const real expterm = std::exp(-(std::max(static_cast<real>(0.0), localSlipRate) + vStar) *
                                  deltaT / prakashLength);

    const real newStrength = ctx.data->regularizedStrength[ctx.ltsFace][ctx.pointIndex] * expterm +
                             faultStrength * (1.0 - expterm);

    ctx.data->regularizedStrength[ctx.ltsFace][ctx.pointIndex] = newStrength;
    return newStrength;
  };
};

class TPApprox {
  public:
  TPApprox(seissol::initializer::parameters::DRParameters* parameters) {};

  static void copyLtsTreeToLocal(FrictionLawData* data,
                                 seissol::initializer::Layer& layerData,
                                 const seissol::initializer::DynamicRupture* const dynRup,
                                 real fullUpdateTime) {}

  SEISSOL_DEVICE static real
      resampleSlipRate(FrictionLawContext& ctx,
                       const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints]) {
    return slipRateMagnitude[ctx.pointIndex];
  };

  SEISSOL_DEVICE static real stateVariableHook(FrictionLawContext& ctx,
                                               real localAccumulatedSlip,
                                               real localDc,
                                               real tpProxyExponent) {
    const real factor = (1.0 + std::fabs(localAccumulatedSlip) / localDc);
    return 1.0 - std::pow(factor, -tpProxyExponent);
  };

  SEISSOL_DEVICE static real strengthHook(FrictionLawContext& ctx,
                                          real strength,
                                          real localSlipRate,
                                          real deltaT,
                                          real vStar,
                                          real prakashLength) {
    return strength;
  };
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_LINEARSLIPWEAKENING_H_
