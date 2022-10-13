#ifndef SEISSOL_GPU_LINEARSLIPWEAKENING_H
#define SEISSOL_GPU_LINEARSLIPWEAKENING_H

#include "DynamicRupture/FrictionLaws/GpuImpl/GpuFrictionSolver.h"

namespace seissol::dr::friction_law::gpu {

/**
 * Abstract Class implementing the general structure of linear slip weakening friction laws.
 * specific implementation is done by overriding and implementing the hook functions (via CRTP).
 */
template <typename Derived>
class LinearSlipWeakeningBase : public GpuFrictionSolver<LinearSlipWeakeningBase<Derived>> {
  public:
  LinearSlipWeakeningBase<Derived>(dr::DRParameters* drParameters)
      : GpuFrictionSolver<LinearSlipWeakeningBase<Derived>>(drParameters){};

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
  void calcSlipRateAndTraction(FaultStresses* faultStressesPtr,
                               TractionResults* tractionResultsPtr,
                               real (*strengthBuffer)[misc::numPaddedPoints],
                               unsigned int timeIndex) {

    auto layerSize{this->currLayerSize};
    auto* initialStressInFaultCS{this->initialStressInFaultCS};
    auto* impAndEta{this->impAndEta};
    auto* slipRateMagnitude{this->slipRateMagnitude};
    auto* slipRate1{this->slipRate1};
    auto* slipRate2{this->slipRate2};
    auto* traction1{this->traction1};
    auto* traction2{this->traction2};
    auto* slip1{this->slip1};
    auto* slip2{this->slip2};
    auto deltaT{this->deltaT[timeIndex]};

    sycl::nd_range rng{{layerSize * misc::numPaddedPoints}, {misc::numPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        auto& faultStresses = faultStressesPtr[ltsFace];
        auto& tractionResults = tractionResultsPtr[ltsFace];
        auto& strength = strengthBuffer[ltsFace];

        // calculate absolute value of stress in Y and Z direction
        const real totalStress1 = initialStressInFaultCS[ltsFace][pointIndex][3] +
                                  faultStresses.traction1[timeIndex][pointIndex];
        const real totalStress2 = initialStressInFaultCS[ltsFace][pointIndex][5] +
                                  faultStresses.traction2[timeIndex][pointIndex];
        const real absoluteShearStress = misc::magnitude(totalStress1, totalStress2);
        // calculate slip rates
        slipRateMagnitude[ltsFace][pointIndex] =
            std::max(static_cast<real>(0.0),
                     (absoluteShearStress - strength[pointIndex]) * impAndEta[ltsFace].invEtaS);
        const auto divisor =
            strength[pointIndex] + impAndEta[ltsFace].etaS * slipRateMagnitude[ltsFace][pointIndex];
        slipRate1[ltsFace][pointIndex] =
            slipRateMagnitude[ltsFace][pointIndex] * totalStress1 / divisor;
        slipRate2[ltsFace][pointIndex] =
            slipRateMagnitude[ltsFace][pointIndex] * totalStress2 / divisor;
        // calculate traction
        tractionResults.traction1[timeIndex][pointIndex] =
            faultStresses.traction1[timeIndex][pointIndex] -
            impAndEta[ltsFace].etaS * slipRate1[ltsFace][pointIndex];
        tractionResults.traction2[timeIndex][pointIndex] =
            faultStresses.traction2[timeIndex][pointIndex] -
            impAndEta[ltsFace].etaS * slipRate2[ltsFace][pointIndex];
        traction1[ltsFace][pointIndex] = tractionResults.traction1[timeIndex][pointIndex];
        traction2[ltsFace][pointIndex] = tractionResults.traction2[timeIndex][pointIndex];
        // update directional slip
        slip1[ltsFace][pointIndex] += slipRate1[ltsFace][pointIndex] * deltaT;
        slip2[ltsFace][pointIndex] += slipRate2[ltsFace][pointIndex] * deltaT;
      });
    });
  }

  /**
   * evaluate friction law: updated mu -> friction law
   * for example see Carsten Uphoff's thesis: Eq. 2.45
   */
  void frictionFunctionHook(real (*stateVariableBuffer)[misc::numPaddedPoints]) {
    auto layerSize{this->currLayerSize};
    auto* mu{this->mu};
    auto* muS{this->muS};
    auto* muD{this->muD};

    sycl::nd_range rng{{layerSize * misc::numPaddedPoints}, {misc::numPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        auto& stateVariable = stateVariableBuffer[ltsFace];
        mu[ltsFace][pointIndex] =
            muS[ltsFace][pointIndex] -
            (muS[ltsFace][pointIndex] - muD[ltsFace][pointIndex]) * stateVariable[pointIndex];
      });
    });
  }

  /*
   * output time when shear stress is equal to the dynamic stress after rupture arrived
   * currently only for linear slip weakening
   */
  void saveDynamicStressOutput() {
    auto layerSize{this->currLayerSize};
    auto fullUpdateTime{this->mFullUpdateTime};
    auto* dynStressTime{this->dynStressTime};
    auto* dynStressTimePending{this->dynStressTimePending};
    auto* accumulatedSlipMagnitude{this->accumulatedSlipMagnitude};
    auto* dC{this->dC};

    sycl::nd_range rng{{layerSize * misc::numPaddedPoints}, {misc::numPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        if (dynStressTimePending[ltsFace][pointIndex] &&
            std::fabs(accumulatedSlipMagnitude[ltsFace][pointIndex]) >= dC[ltsFace][pointIndex]) {
          dynStressTime[ltsFace][pointIndex] = fullUpdateTime;
          dynStressTimePending[ltsFace][pointIndex] = false;
        }
      });
    });
  }

  void preHook(real (*stateVariableBuffer)[misc::numPaddedPoints]) {}
  void postHook(real (*stateVariableBuffer)[misc::numPaddedPoints]) {}

  protected:
  static constexpr real u0 = 10e-14;
  real (*dC)[misc::numPaddedPoints];
  real (*muS)[misc::numPaddedPoints];
  real (*muD)[misc::numPaddedPoints];
  real (*cohesion)[misc::numPaddedPoints];
  real (*forcedRuptureTime)[misc::numPaddedPoints];

  // Note, regularisedStrength belongs only to `BiMaterialFault` but
  // we declare it here to avoid an implicit `this` pointer in the
  // kernel body
  real (*regularisedStrength)[misc::numPaddedPoints];
};

template <class SpecializationT>
class LinearSlipWeakeningLaw
    : public LinearSlipWeakeningBase<LinearSlipWeakeningLaw<SpecializationT>> {
  public:
  LinearSlipWeakeningLaw<SpecializationT>(dr::DRParameters* drParameters)
      : LinearSlipWeakeningBase<LinearSlipWeakeningLaw<SpecializationT>>(drParameters),
        specialization(drParameters){};

  void copySpecificLtsDataTreeToLocal(seissol::initializers::Layer& layerData,
                                      seissol::initializers::DynamicRupture const* const dynRup,
                                      real fullUpdateTime) override {
    auto* concreteLts =
        dynamic_cast<seissol::initializers::LTSLinearSlipWeakening const* const>(dynRup);
    this->dC = layerData.var(concreteLts->dC);
    this->muS = layerData.var(concreteLts->muS);
    this->muD = layerData.var(concreteLts->muD);
    this->cohesion = layerData.var(concreteLts->cohesion);
    this->forcedRuptureTime = layerData.var(concreteLts->forcedRuptureTime);
    this->regularisedStrength = specialization.getRegularisedStrength(layerData, dynRup);
  }

  void calcStrengthHook(FaultStresses* faultStressesPtr,
                        real (*strengthBuffer)[misc::numPaddedPoints],
                        unsigned int timeIndex) {

    auto layerSize{this->currLayerSize};
    auto deltaT{this->deltaT[timeIndex]};
    auto* initialStressInFaultCS{this->initialStressInFaultCS};
    auto* slipRateMagnitude{this->slipRateMagnitude};
    auto* cohesion{this->cohesion};
    auto* mu{this->mu};

    auto* regularisedStrength{this->regularisedStrength};
    const auto vStar{this->drParameters->vStar};
    const auto prakashLength{this->drParameters->prakashLength};

    sycl::nd_range rng{{layerSize * misc::numPaddedPoints}, {misc::numPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        auto& faultStresses = faultStressesPtr[ltsFace];
        auto& strength = strengthBuffer[ltsFace];

        const real totalNormalStress = initialStressInFaultCS[ltsFace][pointIndex][0] +
                                       faultStresses.normalStress[timeIndex][pointIndex];
        strength[pointIndex] =
            -cohesion[ltsFace][pointIndex] -
            mu[ltsFace][pointIndex] * std::min(totalNormalStress, static_cast<real>(0.0));

        strength[pointIndex] = SpecializationT::strengthHook(regularisedStrength[ltsFace],
                                                             strength[pointIndex],
                                                             slipRateMagnitude[ltsFace][pointIndex],
                                                             deltaT,
                                                             vStar,
                                                             prakashLength,
                                                             pointIndex);
      });
    });
  }

  void calcStateVariableHook(real (*stateVariableBuffer)[misc::numPaddedPoints],
                             unsigned int timeIndex) {

    auto layerSize{this->currLayerSize};
    auto* accumulatedSlipMagnitude{this->accumulatedSlipMagnitude};
    auto* slipRateMagnitude{this->slipRateMagnitude};
    auto* forcedRuptureTime{this->forcedRuptureTime};
    auto* dC{this->dC};
    auto* resample{this->resampleMatrix};
    auto deltaT{this->deltaT[timeIndex]};
    const real tn{this->mFullUpdateTime + deltaT};
    const auto t0{this->drParameters->t0};

    sycl::nd_range rng{{layerSize * misc::numPaddedPoints}, {misc::numPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      sycl::accessor<real, 1, sycl::access::mode::read_write, sycl::access::target::local>
          resampledSlipRate(misc::numPaddedPoints, cgh);
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        auto ltsFace = item.get_group().get_group_id(0);
        auto pointIndex = item.get_local_id(0);

        resampledSlipRate[pointIndex] =
            specialization.resampleSlipRate(resample, slipRateMagnitude[ltsFace], pointIndex);

        item.barrier(sycl::access::fence_space::local_space);
        auto& stateVariable = stateVariableBuffer[ltsFace];

        // integrate slip rate to get slip = state variable
        accumulatedSlipMagnitude[ltsFace][pointIndex] += resampledSlipRate[pointIndex] * deltaT;

        // Modif T. Ulrich-> generalisation of tpv16/17 to 30/31
        // Actually slip is already the stateVariable for this FL, but to simplify the next
        // equations we divide it here by the critical distance.
        const real localStateVariable = std::min(
            std::fabs(accumulatedSlipMagnitude[ltsFace][pointIndex]) / dC[ltsFace][pointIndex],
            static_cast<real>(1.0));

        real f2 = 0.0;
        if (t0 == 0) {
          f2 = 1.0 * (tn >= forcedRuptureTime[ltsFace][pointIndex]);
        } else {
          f2 = std::clamp((tn - forcedRuptureTime[ltsFace][pointIndex]) / t0,
                          static_cast<real>(0.0),
                          static_cast<real>(1.0));
        }
        stateVariable[pointIndex] = std::max(localStateVariable, f2);
      });
    });
  }

  protected:
  SpecializationT specialization;
};

class NoSpecialization {
  public:
  NoSpecialization(DRParameters* parameters){};

  using returnArrayType = real (*)[dr::misc::numPaddedPoints];
  returnArrayType
      getRegularisedStrength(seissol::initializers::Layer& layerData,
                             seissol::initializers::DynamicRupture const* const dynRup) {
    // `NoSpecialization` type does not provide `regularisedStrength`
    return nullptr;
  }

  real resampleSlipRate(real const* resampleMatrix,
                        real const (&slipRateMagnitude)[dr::misc::numPaddedPoints],
                        size_t pointIndex) {

    // perform matrix vector multiplication

    constexpr auto dim0 = misc::dimSize<init::resample, 0>();
    constexpr auto dim1 = misc::dimSize<init::resample, 1>();
    static_assert(dim0 == misc::numPaddedPoints);
    static_assert(dim0 >= dim1);

    real result{0.0};
    for (size_t i{0}; i < dim1; ++i) {
      result += resampleMatrix[pointIndex + i * dim0] * slipRateMagnitude[i];
    }
    return result;
  };

  static real strengthHook(real regularisedStrength[dr::misc::numPaddedPoints],
                           real strength,
                           real localSlipRate,
                           real deltaT,
                           real vStar,
                           real prakashLength,
                           size_t pointIndex) {
    return strength;
  };
};

class BiMaterialFault {
  public:
  BiMaterialFault(DRParameters* parameters){};

  auto getRegularisedStrength(seissol::initializers::Layer& layerData,
                              seissol::initializers::DynamicRupture const* const dynRup) {
    auto* concreteLts =
        dynamic_cast<seissol::initializers::LTSLinearSlipWeakeningBimaterial const* const>(dynRup);
    return layerData.var(concreteLts->regularisedStrength);
  }

  real resampleSlipRate([[maybe_unused]] real const* resampleMatrix,
                        real const (&slipRateMagnitude)[dr::misc::numPaddedPoints],
                        size_t pointIndex) {
    return slipRateMagnitude[pointIndex];
  };

  static real strengthHook(real regularisedStrength[dr::misc::numPaddedPoints],
                           real faultStrength,
                           real localSlipRate,
                           real deltaT,
                           real vStar,
                           real prakashLength,
                           size_t pointIndex) {

    assert(regularisedStrength != nullptr && "regularisedStrength is not initialized");

    const real expterm = sycl::exp(-(sycl::max(static_cast<real>(0.0), localSlipRate) + vStar) *
                                   deltaT / prakashLength);

    const real newStrength =
        regularisedStrength[pointIndex] * expterm + faultStrength * (1.0 - expterm);

    regularisedStrength[pointIndex] = newStrength;
    return newStrength;
  };
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPU_LINEARSLIPWEAKENING_H
