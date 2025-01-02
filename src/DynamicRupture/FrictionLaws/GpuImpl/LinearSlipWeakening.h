#ifndef SEISSOL_GPU_LINEARSLIPWEAKENING_H
#define SEISSOL_GPU_LINEARSLIPWEAKENING_H

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"

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

  void updateFrictionAndSlip(unsigned timeIndex) {
    // computes fault strength, which is the critical value whether active slip exists.
    static_cast<Derived*>(this)->calcStrengthHook(
        this->faultStresses, this->strengthBuffer, timeIndex);
    // computes resulting slip rates, traction and slip dependent on current friction
    // coefficient and strength
    this->calcSlipRateAndTraction(
        this->faultStresses, this->tractionResults, this->strengthBuffer, timeIndex);
    static_cast<Derived*>(this)->calcStateVariableHook(this->stateVariableBuffer, timeIndex);
    this->frictionFunctionHook(this->stateVariableBuffer);
  }

  /**
   *  compute the slip rate and the traction from the fault strength and fault stresses
   *  also updates the directional slip1 and slip2
   */
  void calcSlipRateAndTraction(FaultStresses* devFaultStresses,
                               TractionResults* devTractionResults,
                               real (*devStrengthBuffer)[misc::NumPaddedPoints],
                               unsigned int timeIndex) {

    auto* devInitialStressInFaultCS{this->initialStressInFaultCS};
    auto* devImpAndEta{this->impAndEta};
    auto* devSlipRateMagnitude{this->slipRateMagnitude};
    auto* devSlipRate1{this->slipRate1};
    auto* devSlipRate2{this->slipRate2};
    auto* devTraction1{this->traction1};
    auto* devTraction2{this->traction2};
    auto* devSlip1{this->slip1};
    auto* devSlip2{this->slip2};
    auto deltaT{this->deltaT[timeIndex]};

    sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints}, {misc::NumPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        auto& faultStresses = devFaultStresses[ltsFace];
        auto& tractionResults = devTractionResults[ltsFace];
        auto& strength = devStrengthBuffer[ltsFace];

        // calculate absolute value of stress in Y and Z direction
        const real totalStress1 = devInitialStressInFaultCS[ltsFace][pointIndex][3] +
                                  faultStresses.traction1[timeIndex][pointIndex];
        const real totalStress2 = devInitialStressInFaultCS[ltsFace][pointIndex][5] +
                                  faultStresses.traction2[timeIndex][pointIndex];
        const real absoluteShearStress = misc::magnitude(totalStress1, totalStress2);
        // calculate slip rates
        devSlipRateMagnitude[ltsFace][pointIndex] =
            sycl::max(static_cast<real>(0.0),
                      (absoluteShearStress - strength[pointIndex]) * devImpAndEta[ltsFace].invEtaS);
        const auto divisor = strength[pointIndex] +
                             devImpAndEta[ltsFace].etaS * devSlipRateMagnitude[ltsFace][pointIndex];
        devSlipRate1[ltsFace][pointIndex] =
            devSlipRateMagnitude[ltsFace][pointIndex] * totalStress1 / divisor;
        devSlipRate2[ltsFace][pointIndex] =
            devSlipRateMagnitude[ltsFace][pointIndex] * totalStress2 / divisor;
        // calculate traction
        tractionResults.traction1[timeIndex][pointIndex] =
            faultStresses.traction1[timeIndex][pointIndex] -
            devImpAndEta[ltsFace].etaS * devSlipRate1[ltsFace][pointIndex];
        tractionResults.traction2[timeIndex][pointIndex] =
            faultStresses.traction2[timeIndex][pointIndex] -
            devImpAndEta[ltsFace].etaS * devSlipRate2[ltsFace][pointIndex];
        devTraction1[ltsFace][pointIndex] = tractionResults.traction1[timeIndex][pointIndex];
        devTraction2[ltsFace][pointIndex] = tractionResults.traction2[timeIndex][pointIndex];
        // update directional slip
        devSlip1[ltsFace][pointIndex] += devSlipRate1[ltsFace][pointIndex] * deltaT;
        devSlip2[ltsFace][pointIndex] += devSlipRate2[ltsFace][pointIndex] * deltaT;
      });
    });
  }

  /**
   * evaluate friction law: updated mu -> friction law
   * for example see Carsten Uphoff's thesis: Eq. 2.45
   */
  void frictionFunctionHook(real (*stateVariableBuffer)[misc::NumPaddedPoints]) {
    auto* devMu{this->mu};
    auto* devMuS{this->muS};
    auto* devMuD{this->muD};
    auto* devSlipRateMagnitude{this->slipRateMagnitude};
    auto* devPeakSlipRate{this->peakSlipRate};
    auto devHealingThreshold{this->drParameters->healingThreshold};

    sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints}, {misc::NumPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        auto& stateVariable = stateVariableBuffer[ltsFace];
        devMu[ltsFace][pointIndex] =
            devMuS[ltsFace][pointIndex] -
            (devMuS[ltsFace][pointIndex] - devMuD[ltsFace][pointIndex]) * stateVariable[pointIndex];
        // instantaneous healing
        if ((devPeakSlipRate[ltsFace][pointIndex] > devHealingThreshold) &&
            (devSlipRateMagnitude[ltsFace][pointIndex] < devHealingThreshold)) {
          devMu[ltsFace][pointIndex] = devMuS[ltsFace][pointIndex];
          stateVariable[pointIndex] = 0.0;
        }
      });
    });
  }

  /*
   * output time when shear stress is equal to the dynamic stress after rupture arrived
   * currently only for linear slip weakening
   */
  void saveDynamicStressOutput() {
    auto fullUpdateTime{this->mFullUpdateTime};
    auto* devDynStressTime{this->dynStressTime};
    auto* devDynStressTimePending{this->dynStressTimePending};
    auto* devAccumulatedSlipMagnitude{this->accumulatedSlipMagnitude};
    auto* devDC{this->dC};

    sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints}, {misc::NumPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        if (devDynStressTimePending[ltsFace][pointIndex] &&
            sycl::fabs(devAccumulatedSlipMagnitude[ltsFace][pointIndex]) >=
                devDC[ltsFace][pointIndex]) {
          devDynStressTime[ltsFace][pointIndex] = fullUpdateTime;
          devDynStressTimePending[ltsFace][pointIndex] = false;
        }
      });
    });
  }

  void preHook(real (*stateVariableBuffer)[misc::NumPaddedPoints]) {}
  void postHook(real (*stateVariableBuffer)[misc::NumPaddedPoints]) {}

  protected:
  static constexpr real U0 = 10e-14;
  real (*dC)[misc::NumPaddedPoints];
  real (*muS)[misc::NumPaddedPoints];
  real (*muD)[misc::NumPaddedPoints];
  real (*cohesion)[misc::NumPaddedPoints];
  real (*forcedRuptureTime)[misc::NumPaddedPoints];
};

template <class SpecializationT>
class LinearSlipWeakeningLaw
    : public LinearSlipWeakeningBase<LinearSlipWeakeningLaw<SpecializationT>> {
  public:
  LinearSlipWeakeningLaw<SpecializationT>(
      seissol::initializer::parameters::DRParameters* drParameters)
      : LinearSlipWeakeningBase<LinearSlipWeakeningLaw<SpecializationT>>(drParameters),
        specialization(drParameters){};

  void copySpecificLtsDataTreeToLocal(seissol::initializer::Layer& layerData,
                                      const seissol::initializer::DynamicRupture* const dynRup,
                                      real fullUpdateTime) override {
    auto* concreteLts = dynamic_cast<const seissol::initializer::LTSLinearSlipWeakening*>(dynRup);
    this->dC = layerData.var(concreteLts->dC, seissol::initializer::AllocationPlace::Device);
    this->muS = layerData.var(concreteLts->muS, seissol::initializer::AllocationPlace::Device);
    this->muD = layerData.var(concreteLts->muD, seissol::initializer::AllocationPlace::Device);
    this->cohesion =
        layerData.var(concreteLts->cohesion, seissol::initializer::AllocationPlace::Device);
    this->forcedRuptureTime = layerData.var(concreteLts->forcedRuptureTime,
                                            seissol::initializer::AllocationPlace::Device);
    this->specialization.copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  }

  void calcStrengthHook(FaultStresses* devFaultStresses,
                        real (*devStrengthBuffer)[misc::NumPaddedPoints],
                        unsigned int timeIndex) {

    auto deltaT{this->deltaT[timeIndex]};
    auto* devInitialStressInFaultCS{this->initialStressInFaultCS};
    auto* devSlipRateMagnitude{this->slipRateMagnitude};
    auto* devCohesion{this->cohesion};
    auto* devMu{this->mu};

    const auto vStar{this->drParameters->vStar};
    const auto prakashLength{this->drParameters->prakashLength};
    auto currentLayerDetails = specialization.getCurrentLayerDetails();

    sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints}, {misc::NumPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        auto& faultStresses = devFaultStresses[ltsFace];
        auto& strength = devStrengthBuffer[ltsFace];

        const real totalNormalStress = devInitialStressInFaultCS[ltsFace][pointIndex][0] +
                                       faultStresses.normalStress[timeIndex][pointIndex];
        strength[pointIndex] =
            -devCohesion[ltsFace][pointIndex] -
            devMu[ltsFace][pointIndex] * sycl::min(totalNormalStress, static_cast<real>(0.0));

        strength[pointIndex] =
            SpecializationT::strengthHook(currentLayerDetails,
                                          strength[pointIndex],
                                          devSlipRateMagnitude[ltsFace][pointIndex],
                                          deltaT,
                                          vStar,
                                          prakashLength,
                                          ltsFace,
                                          pointIndex);
      });
    });
  }

  void calcStateVariableHook(real (*devStateVariableBuffer)[misc::NumPaddedPoints],
                             unsigned int timeIndex) {

    auto* devAccumulatedSlipMagnitude{this->accumulatedSlipMagnitude};
    auto* devSlipRateMagnitude{this->slipRateMagnitude};
    auto* devForcedRuptureTime{this->forcedRuptureTime};
    auto* devDC{this->dC};
    auto* devResample{this->resampleMatrix};
    auto deltaT{this->deltaT[timeIndex]};
    const real tn{this->mFullUpdateTime + deltaT};
    const auto t0{this->drParameters->t0};
    const auto tpProxyExponent{this->drParameters->tpProxyExponent};

    sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints}, {misc::NumPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        auto ltsFace = item.get_group().get_group_id(0);
        auto pointIndex = item.get_local_id(0);

        const real resampledSlipRate = SpecializationT::resampleSlipRate(
            devResample, devSlipRateMagnitude[ltsFace], pointIndex);

        // integrate slip rate to get slip = state variable
        devAccumulatedSlipMagnitude[ltsFace][pointIndex] += resampledSlipRate * deltaT;

        // Actually slip is already the stateVariable for this FL, but to simplify the next
        // equations we divide it here by the critical distance.
        const real localStateVariable =
            SpecializationT::stateVariableHook(devAccumulatedSlipMagnitude[ltsFace][pointIndex],
                                               devDC[ltsFace][pointIndex],
                                               tpProxyExponent,
                                               ltsFace,
                                               pointIndex);

        real f2 = 0.0;
        if (t0 == 0) {
          f2 = 1.0 * (tn >= devForcedRuptureTime[ltsFace][pointIndex]);
        } else {
          f2 = sycl::clamp((tn - devForcedRuptureTime[ltsFace][pointIndex]) / t0,
                           static_cast<real>(0.0),
                           static_cast<real>(1.0));
        }
        devStateVariableBuffer[ltsFace][pointIndex] = sycl::max(localStateVariable, f2);
      });
    });
  }

  protected:
  SpecializationT specialization;
};

class NoSpecialization {
  public:
  NoSpecialization(seissol::initializer::parameters::DRParameters* parameters) {};

  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime) {}

  static real resampleSlipRate(const real* resampleMatrix,
                               const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints],
                               size_t pointIndex) {

    // perform matrix vector multiplication

    constexpr auto Dim0 = misc::dimSize<init::resample, 0>();
    constexpr auto Dim1 = misc::dimSize<init::resample, 1>();
    static_assert(Dim0 == misc::NumPaddedPoints);
    static_assert(Dim0 >= Dim1);

    real result{0.0};
    for (size_t i{0}; i < Dim1; ++i) {
      result += resampleMatrix[pointIndex + i * Dim0] * slipRateMagnitude[i];
    }
    return result;
  };

  struct Details {};
  Details getCurrentLayerDetails() { return Details{}; }

  static real stateVariableHook(real localAccumulatedSlip,
                                real localDc,
                                real tpProxyExponent,
                                size_t ltsFace,
                                size_t pointIndex) {
    return sycl::min(sycl::fabs(localAccumulatedSlip) / localDc, static_cast<real>(1.0));
  };

  static real strengthHook(Details details,
                           real strength,
                           real localSlipRate,
                           real deltaT,
                           real vStar,
                           real prakashLength,
                           size_t ltsFace,
                           size_t pointIndex) {
    return strength;
  };
};

class BiMaterialFault {
  public:
  BiMaterialFault(seissol::initializer::parameters::DRParameters* parameters) {};

  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime) {
    auto* concreteLts =
        dynamic_cast<const seissol::initializer::LTSLinearSlipWeakeningBimaterial*>(dynRup);
    this->regularisedStrength = layerData.var(concreteLts->regularisedStrength,
                                              seissol::initializer::AllocationPlace::Device);
  }

  static real resampleSlipRate([[maybe_unused]] const real* resampleMatrix,
                               const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints],
                               size_t pointIndex) {
    return slipRateMagnitude[pointIndex];
  };

  struct Details {
    real (*regularisedStrength)[misc::NumPaddedPoints];
  };

  Details getCurrentLayerDetails() {
    Details details{this->regularisedStrength};
    return details;
  }

  static real stateVariableHook(real localAccumulatedSlip,
                                real localDc,
                                real tpProxyExponent,
                                size_t ltsFace,
                                size_t pointIndex) {
    return sycl::min(sycl::fabs(localAccumulatedSlip) / localDc, static_cast<real>(1.0));
  };

  static real strengthHook(Details details,
                           real faultStrength,
                           real localSlipRate,
                           real deltaT,
                           real vStar,
                           real prakashLength,
                           size_t ltsFace,
                           size_t pointIndex) {

    auto* regularisedStrength = details.regularisedStrength[ltsFace];

    const real expterm = sycl::exp(-(sycl::max(static_cast<real>(0.0), localSlipRate) + vStar) *
                                   deltaT / prakashLength);

    const real newStrength =
        regularisedStrength[pointIndex] * expterm + faultStrength * (1.0 - expterm);

    regularisedStrength[pointIndex] = newStrength;
    return newStrength;
  };

  private:
  real (*regularisedStrength)[misc::NumPaddedPoints];
};

class TPApprox {
  public:
  TPApprox(seissol::initializer::parameters::DRParameters* parameters) {};

  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime) {}

  static real resampleSlipRate([[maybe_unused]] const real* resampleMatrix,
                               const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints],
                               size_t pointIndex) {
    return slipRateMagnitude[pointIndex];
  };

  struct Details {
    real (*regularisedStrength)[misc::NumPaddedPoints];
  };

  Details getCurrentLayerDetails() { return Details{}; }

  static real stateVariableHook(real localAccumulatedSlip,
                                real localDc,
                                real tpProxyExponent,
                                size_t ltsFace,
                                size_t pointIndex) {
    const real factor = (1.0 + sycl::fabs(localAccumulatedSlip) / localDc);
    return 1.0 - sycl::pow(factor, -tpProxyExponent);
  };

  static real strengthHook(Details details,
                           real strength,
                           real localSlipRate,
                           real deltaT,
                           real vStar,
                           real prakashLength,
                           size_t ltsFace,
                           size_t pointIndex) {
    return strength;
  };
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPU_LINEARSLIPWEAKENING_H
