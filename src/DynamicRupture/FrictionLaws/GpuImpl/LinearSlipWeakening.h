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

  static void updateFrictionAndSlip(FrictionLawContext& ctx, unsigned timeIndex) {
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
  static void calcSlipRateAndTraction(FrictionLawContext& ctx, unsigned int timeIndex) {

    auto* devInitialStressInFaultCS{ctx.data->initialStressInFaultCS};
    auto* devImpAndEta{ctx.data->impAndEta};
    auto* devSlipRateMagnitude{ctx.data->slipRateMagnitude};
    auto* devSlipRate1{ctx.data->slipRate1};
    auto* devSlipRate2{ctx.data->slipRate2};
    auto* devTraction1{ctx.data->traction1};
    auto* devTraction2{ctx.data->traction2};
    auto* devSlip1{ctx.data->slip1};
    auto* devSlip2{ctx.data->slip2};
    auto deltaT{ctx.data->deltaT[timeIndex]};

    auto& faultStresses = ctx.faultStresses;
    auto& tractionResults = ctx.tractionResults;
    auto& strength = ctx.strengthBuffer;

    // calculate absolute value of stress in Y and Z direction
    const real totalStress1 = devInitialStressInFaultCS[ctx.ltsFace][ctx.pointIndex][3] +
                              faultStresses.traction1[timeIndex];
    const real totalStress2 = devInitialStressInFaultCS[ctx.ltsFace][ctx.pointIndex][5] +
                              faultStresses.traction2[timeIndex];
    const real absoluteShearStress = misc::magnitude(totalStress1, totalStress2);
    // calculate slip rates
    devSlipRateMagnitude[ctx.ltsFace][ctx.pointIndex] =
        sycl::max(static_cast<real>(0.0),
                  (absoluteShearStress - strength) * devImpAndEta[ctx.ltsFace].invEtaS);
    const auto divisor = strength + devImpAndEta[ctx.ltsFace].etaS *
                                        devSlipRateMagnitude[ctx.ltsFace][ctx.pointIndex];
    devSlipRate1[ctx.ltsFace][ctx.pointIndex] =
        devSlipRateMagnitude[ctx.ltsFace][ctx.pointIndex] * totalStress1 / divisor;
    devSlipRate2[ctx.ltsFace][ctx.pointIndex] =
        devSlipRateMagnitude[ctx.ltsFace][ctx.pointIndex] * totalStress2 / divisor;
    // calculate traction
    tractionResults.traction1[timeIndex] =
        faultStresses.traction1[timeIndex] -
        devImpAndEta[ctx.ltsFace].etaS * devSlipRate1[ctx.ltsFace][ctx.pointIndex];
    tractionResults.traction2[timeIndex] =
        faultStresses.traction2[timeIndex] -
        devImpAndEta[ctx.ltsFace].etaS * devSlipRate2[ctx.ltsFace][ctx.pointIndex];
    devTraction1[ctx.ltsFace][ctx.pointIndex] = tractionResults.traction1[timeIndex];
    devTraction2[ctx.ltsFace][ctx.pointIndex] = tractionResults.traction2[timeIndex];
    // update directional slip
    devSlip1[ctx.ltsFace][ctx.pointIndex] += devSlipRate1[ctx.ltsFace][ctx.pointIndex] * deltaT;
    devSlip2[ctx.ltsFace][ctx.pointIndex] += devSlipRate2[ctx.ltsFace][ctx.pointIndex] * deltaT;
  }

  /**
   * evaluate friction law: updated mu -> friction law
   * for example see Carsten Uphoff's thesis: Eq. 2.45
   */
  static void frictionFunctionHook(FrictionLawContext& ctx) {
    auto* devMu{ctx.data->mu};
    auto* devMuS{ctx.data->muS};
    auto* devMuD{ctx.data->muD};
    auto* devSlipRateMagnitude{ctx.data->slipRateMagnitude};
    auto* devPeakSlipRate{ctx.data->peakSlipRate};
    auto devHealingThreshold{ctx.data->drParameters.healingThreshold};

    auto& stateVariable = ctx.stateVariableBuffer;
    devMu[ctx.ltsFace][ctx.pointIndex] =
        devMuS[ctx.ltsFace][ctx.pointIndex] -
        (devMuS[ctx.ltsFace][ctx.pointIndex] - devMuD[ctx.ltsFace][ctx.pointIndex]) * stateVariable;
    // instantaneous healing
    if ((devPeakSlipRate[ctx.ltsFace][ctx.pointIndex] > devHealingThreshold) &&
        (devSlipRateMagnitude[ctx.ltsFace][ctx.pointIndex] < devHealingThreshold)) {
      devMu[ctx.ltsFace][ctx.pointIndex] = devMuS[ctx.ltsFace][ctx.pointIndex];
      stateVariable = 0.0;
    }
  }

  /*
   * output time when shear stress is equal to the dynamic stress after rupture arrived
   * currently only for linear slip weakening
   */
  static void saveDynamicStressOutput(FrictionLawContext& ctx) {
    auto fullUpdateTime{ctx.data->mFullUpdateTime};
    auto* devDynStressTime{ctx.data->dynStressTime};
    auto* devDynStressTimePending{ctx.data->dynStressTimePending};
    auto* devAccumulatedSlipMagnitude{ctx.data->accumulatedSlipMagnitude};
    auto* devDC{ctx.data->dC};

    if (devDynStressTimePending[ctx.ltsFace][ctx.pointIndex] &&
        sycl::fabs(devAccumulatedSlipMagnitude[ctx.ltsFace][ctx.pointIndex]) >=
            devDC[ctx.ltsFace][ctx.pointIndex]) {
      devDynStressTime[ctx.ltsFace][ctx.pointIndex] = fullUpdateTime;
      devDynStressTimePending[ctx.ltsFace][ctx.pointIndex] = false;
    }
  }

  static void preHook(FrictionLawContext& ctx) {}
  static void postHook(FrictionLawContext& ctx) {}

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

  static void calcStrengthHook(FrictionLawContext& ctx, unsigned int timeIndex) {

    auto deltaT{ctx.data->deltaT[timeIndex]};
    auto* devInitialStressInFaultCS{ctx.data->initialStressInFaultCS};
    auto* devSlipRateMagnitude{ctx.data->slipRateMagnitude};
    auto* devCohesion{ctx.data->cohesion};
    auto* devMu{ctx.data->mu};

    const auto vStar{ctx.data->drParameters.vStar};
    const auto prakashLength{ctx.data->drParameters.prakashLength};

    auto& faultStresses = ctx.faultStresses;
    auto& strength = ctx.strengthBuffer;

    const real totalNormalStress = devInitialStressInFaultCS[ctx.ltsFace][ctx.pointIndex][0] +
                                   faultStresses.normalStress[timeIndex];
    strength =
        -devCohesion[ctx.ltsFace][ctx.pointIndex] -
        devMu[ctx.ltsFace][ctx.pointIndex] * sycl::min(totalNormalStress, static_cast<real>(0.0));

    strength = SpecializationT::strengthHook(ctx,
                                             strength,
                                             devSlipRateMagnitude[ctx.ltsFace][ctx.pointIndex],
                                             deltaT,
                                             vStar,
                                             prakashLength);
  }

  static void calcStateVariableHook(FrictionLawContext& ctx, unsigned int timeIndex) {

    auto* devAccumulatedSlipMagnitude{ctx.data->accumulatedSlipMagnitude};
    auto* devSlipRateMagnitude{ctx.data->slipRateMagnitude};
    auto* devForcedRuptureTime{ctx.data->forcedRuptureTime};
    auto* devDC{ctx.data->dC};
    auto* devResample{ctx.resampleMatrix};
    auto deltaT{ctx.data->deltaT[timeIndex]};
    const real tn{ctx.data->mFullUpdateTime + deltaT};
    const auto t0{ctx.data->drParameters.t0};
    const auto tpProxyExponent{ctx.data->drParameters.tpProxyExponent};

    const real resampledSlipRate =
        SpecializationT::resampleSlipRate(ctx, devSlipRateMagnitude[ctx.ltsFace]);

    // integrate slip rate to get slip = state variable
    devAccumulatedSlipMagnitude[ctx.ltsFace][ctx.pointIndex] += resampledSlipRate * deltaT;

    // Actually slip is already the stateVariable for this FL, but to simplify the next
    // equations we divide it here by the critical distance.
    const real localStateVariable =
        SpecializationT::stateVariableHook(ctx,
                                           devAccumulatedSlipMagnitude[ctx.ltsFace][ctx.pointIndex],
                                           devDC[ctx.ltsFace][ctx.pointIndex],
                                           tpProxyExponent);

    real f2 = 0.0;
    if (t0 == 0) {
      f2 = 1.0 * static_cast<double>(tn >= devForcedRuptureTime[ctx.ltsFace][ctx.pointIndex]);
    } else {
      f2 = sycl::clamp((tn - devForcedRuptureTime[ctx.ltsFace][ctx.pointIndex]) / t0,
                       static_cast<real>(0.0),
                       static_cast<real>(1.0));
    }
    ctx.stateVariableBuffer = sycl::max(localStateVariable, f2);
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

  static real resampleSlipRate(FrictionLawContext& ctx,
                               const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints]) {

    // perform matrix vector multiplication

    constexpr auto Dim0 = misc::dimSize<init::resample, 0>();
    constexpr auto Dim1 = misc::dimSize<init::resample, 1>();
    static_assert(Dim0 == misc::NumPaddedPoints);
    static_assert(Dim0 >= Dim1);

    real result{0.0};
    for (size_t i{0}; i < Dim1; ++i) {
      result += ctx.resampleMatrix[ctx.pointIndex + i * Dim0] * slipRateMagnitude[i];
    }
    return result;
  };

  static real stateVariableHook(FrictionLawContext& ctx,
                                real localAccumulatedSlip,
                                real localDc,
                                real tpProxyExponent) {
    return sycl::min(sycl::fabs(localAccumulatedSlip) / localDc, static_cast<real>(1.0));
  };

  static real strengthHook(FrictionLawContext& ctx,
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
    data->regularisedStrength = layerData.var(concreteLts->regularisedStrength,
                                              seissol::initializer::AllocationPlace::Device);
  }

  static real resampleSlipRate(FrictionLawContext& ctx,
                               const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints]) {
    return slipRateMagnitude[ctx.pointIndex];
  };

  static real stateVariableHook(FrictionLawContext& ctx,
                                real localAccumulatedSlip,
                                real localDc,
                                real tpProxyExponent) {
    return sycl::min(sycl::fabs(localAccumulatedSlip) / localDc, static_cast<real>(1.0));
  };

  static real strengthHook(FrictionLawContext& ctx,
                           real faultStrength,
                           real localSlipRate,
                           real deltaT,
                           real vStar,
                           real prakashLength) {

    auto* regularisedStrength = ctx.data->regularisedStrength[ctx.ltsFace];

    const real expterm = sycl::exp(-(sycl::max(static_cast<real>(0.0), localSlipRate) + vStar) *
                                   deltaT / prakashLength);

    const real newStrength =
        regularisedStrength[ctx.pointIndex] * expterm + faultStrength * (1.0 - expterm);

    regularisedStrength[ctx.pointIndex] = newStrength;
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

  static real resampleSlipRate(FrictionLawContext& ctx,
                               const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints]) {
    return slipRateMagnitude[ctx.pointIndex];
  };

  static real stateVariableHook(FrictionLawContext& ctx,
                                real localAccumulatedSlip,
                                real localDc,
                                real tpProxyExponent) {
    const real factor = (1.0 + sycl::fabs(localAccumulatedSlip) / localDc);
    return 1.0 - sycl::pow(factor, -tpProxyExponent);
  };

  static real strengthHook(FrictionLawContext& ctx,
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
