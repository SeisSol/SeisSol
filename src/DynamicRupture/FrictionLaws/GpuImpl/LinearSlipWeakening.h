#ifndef SEISSOL_GPU_LINEARSLIPWEAKENING_H
#define SEISSOL_GPU_LINEARSLIPWEAKENING_H

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"

namespace seissol::dr::friction_law::gpu {

/**
 * Abstract Class implementing the general structure of linear slip weakening friction laws.
 * specific implementation is done by overriding and implementing the hook functions (via CRTP).
 */
template <typename Config, typename Derived>
class LinearSlipWeakeningBase
    : public BaseFrictionSolver<Config, LinearSlipWeakeningBase<Config, Derived>> {
  public:
  using RealT = typename Config::RealT;
  LinearSlipWeakeningBase<Config, Derived>(dr::DRParameters* drParameters)
      : BaseFrictionSolver<Config, LinearSlipWeakeningBase<Config, Derived>>(drParameters){};

  void allocateAuxiliaryMemory() override {
    FrictionSolverDetails<Config>::allocateAuxiliaryMemory();
  }

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
  void calcSlipRateAndTraction(FaultStresses<Config>* devFaultStresses,
                               TractionResults<Config>* devTractionResults,
                               RealT (*devStrengthBuffer)[misc::numPaddedPoints<Config>],
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

    sycl::nd_range rng{{this->currLayerSize * misc::numPaddedPoints<Config>},
                       {misc::numPaddedPoints<Config>}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        auto& faultStresses = devFaultStresses[ltsFace];
        auto& tractionResults = devTractionResults[ltsFace];
        auto& strength = devStrengthBuffer[ltsFace];

        // calculate absolute value of stress in Y and Z direction
        const RealT totalStress1 = devInitialStressInFaultCS[ltsFace][pointIndex][3] +
                                   faultStresses.traction1[timeIndex][pointIndex];
        const RealT totalStress2 = devInitialStressInFaultCS[ltsFace][pointIndex][5] +
                                   faultStresses.traction2[timeIndex][pointIndex];
        const RealT absoluteShearStress = misc::magnitude(totalStress1, totalStress2);
        // calculate slip rates
        devSlipRateMagnitude[ltsFace][pointIndex] =
            std::max(static_cast<RealT>(0.0),
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
  void frictionFunctionHook(RealT (*stateVariableBuffer)[misc::numPaddedPoints<Config>]) {
    auto* devMu{this->mu};
    auto* devMuS{this->muS};
    auto* devMuD{this->muD};

    sycl::nd_range rng{{this->currLayerSize * misc::numPaddedPoints<Config>},
                       {misc::numPaddedPoints<Config>}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        auto& stateVariable = stateVariableBuffer[ltsFace];
        devMu[ltsFace][pointIndex] =
            devMuS[ltsFace][pointIndex] -
            (devMuS[ltsFace][pointIndex] - devMuD[ltsFace][pointIndex]) * stateVariable[pointIndex];
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

    sycl::nd_range rng{{this->currLayerSize * misc::numPaddedPoints<Config>},
                       {misc::numPaddedPoints<Config>}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        if (devDynStressTimePending[ltsFace][pointIndex] &&
            std::fabs(devAccumulatedSlipMagnitude[ltsFace][pointIndex]) >=
                devDC[ltsFace][pointIndex]) {
          devDynStressTime[ltsFace][pointIndex] = fullUpdateTime;
          devDynStressTimePending[ltsFace][pointIndex] = false;
        }
      });
    });
  }

  void preHook(RealT (*stateVariableBuffer)[misc::numPaddedPoints<Config>]) {}
  void postHook(RealT (*stateVariableBuffer)[misc::numPaddedPoints<Config>]) {}

  protected:
  static constexpr RealT u0 = 10e-14;
  RealT (*dC)[misc::numPaddedPoints<Config>];
  RealT (*muS)[misc::numPaddedPoints<Config>];
  RealT (*muD)[misc::numPaddedPoints<Config>];
  RealT (*cohesion)[misc::numPaddedPoints<Config>];
  RealT (*forcedRuptureTime)[misc::numPaddedPoints<Config>];
};

template <typename Config, class SpecializationT>
class LinearSlipWeakeningLaw
    : public LinearSlipWeakeningBase<Config, LinearSlipWeakeningLaw<Config, SpecializationT>> {
  public:
  using RealT = typename Config::RealT;
  LinearSlipWeakeningLaw<SpecializationT>(dr::DRParameters* drParameters)
      : LinearSlipWeakeningBase<Config, LinearSlipWeakeningLaw<Config, SpecializationT>>(
            drParameters),
        specialization(drParameters){};

  void copySpecificLtsDataTreeToLocal(
      seissol::initializers::Layer& layerData,
      seissol::initializers::DynamicRupture<Config> const* const dynRup,
      RealT fullUpdateTime) override {
    auto* concreteLts =
        dynamic_cast<seissol::initializers::LTSLinearSlipWeakening<Config> const* const>(dynRup);
    this->dC = layerData.var(concreteLts->dC);
    this->muS = layerData.var(concreteLts->muS);
    this->muD = layerData.var(concreteLts->muD);
    this->cohesion = layerData.var(concreteLts->cohesion);
    this->forcedRuptureTime = layerData.var(concreteLts->forcedRuptureTime);
    this->specialization.copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  }

  void calcStrengthHook(FaultStresses<Config>* devFaultStresses,
                        RealT (*devStrengthBuffer)[misc::numPaddedPoints<Config>],
                        unsigned int timeIndex) {

    auto deltaT{this->deltaT[timeIndex]};
    auto* devInitialStressInFaultCS{this->initialStressInFaultCS};
    auto* devSlipRateMagnitude{this->slipRateMagnitude};
    auto* devCohesion{this->cohesion};
    auto* devMu{this->mu};

    const auto vStar{this->drParameters->vStar};
    const auto prakashLength{this->drParameters->prakashLength};
    auto currentLayerDetails = specialization.getCurrentLayerDetails();

    sycl::nd_range rng{{this->currLayerSize * misc::numPaddedPoints<Config>},
                       {misc::numPaddedPoints<Config>}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        auto& faultStresses = devFaultStresses[ltsFace];
        auto& strength = devStrengthBuffer[ltsFace];

        const RealT totalNormalStress = devInitialStressInFaultCS[ltsFace][pointIndex][0] +
                                        faultStresses.normalStress[timeIndex][pointIndex];
        strength[pointIndex] =
            -devCohesion[ltsFace][pointIndex] -
            devMu[ltsFace][pointIndex] * std::min(totalNormalStress, static_cast<RealT>(0.0));

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

  void calcStateVariableHook(RealT (*devStateVariableBuffer)[misc::numPaddedPoints<Config>],
                             unsigned int timeIndex) {

    auto* devAccumulatedSlipMagnitude{this->accumulatedSlipMagnitude};
    auto* devSlipRateMagnitude{this->slipRateMagnitude};
    auto* devForcedRuptureTime{this->forcedRuptureTime};
    auto* devDC{this->dC};
    auto* devResample{this->resampleMatrix};
    auto deltaT{this->deltaT[timeIndex]};
    const RealT tn{this->mFullUpdateTime + deltaT};
    const auto t0{this->drParameters->t0};

    sycl::nd_range rng{{this->currLayerSize * misc::numPaddedPoints<Config>},
                       {misc::numPaddedPoints<Config>}};
    this->queue.submit([&](sycl::handler& cgh) {
      sycl::accessor<RealT, 1, sycl::access::mode::read_write, sycl::access::target::local>
          resampledSlipRate(misc::numPaddedPoints, cgh);
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        auto ltsFace = item.get_group().get_group_id(0);
        auto pointIndex = item.get_local_id(0);

        resampledSlipRate[pointIndex] = SpecializationT::resampleSlipRate(
            devResample, devSlipRateMagnitude[ltsFace], pointIndex);

        item.barrier(sycl::access::fence_space::local_space);
        auto& stateVariable = devStateVariableBuffer[ltsFace];

        // integrate slip rate to get slip = state variable
        devAccumulatedSlipMagnitude[ltsFace][pointIndex] += resampledSlipRate[pointIndex] * deltaT;

        // Modif T. Ulrich-> generalisation of tpv16/17 to 30/31
        // Actually slip is already the stateVariable for this FL, but to simplify the next
        // equations we divide it here by the critical distance.
        const RealT localStateVariable =
            std::min(std::fabs(devAccumulatedSlipMagnitude[ltsFace][pointIndex]) /
                         devDC[ltsFace][pointIndex],
                     static_cast<RealT>(1.0));

        RealT f2 = 0.0;
        if (t0 == 0) {
          f2 = 1.0 * (tn >= devForcedRuptureTime[ltsFace][pointIndex]);
        } else {
          f2 = std::clamp((tn - devForcedRuptureTime[ltsFace][pointIndex]) / t0,
                          static_cast<RealT>(0.0),
                          static_cast<RealT>(1.0));
        }
        stateVariable[pointIndex] = std::max(localStateVariable, f2);
      });
    });
  }

  protected:
  SpecializationT specialization;
};

template <typename Config>
class NoSpecialization {
  public:
  using RealT = typename Config::RealT;
  NoSpecialization(DRParameters* parameters){};

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture<Config> const* const dynRup,
                          RealT fullUpdateTime) {}

  static RealT resampleSlipRate(RealT const* resampleMatrix,
                                RealT const (&slipRateMagnitude)[dr::misc::numPaddedPoints<Config>],
                                size_t pointIndex) {

    // perform matrix vector multiplication

    constexpr auto dim0 = misc::dimSize<Yateto<Config>::Init::resample, 0>();
    constexpr auto dim1 = misc::dimSize<Yateto<Config>::Init::resample, 1>();
    static_assert(dim0 == misc::numPaddedPoints<Config>);
    static_assert(dim0 >= dim1);

    RealT result{0.0};
    for (size_t i{0}; i < dim1; ++i) {
      result += resampleMatrix[pointIndex + i * dim0] * slipRateMagnitude[i];
    }
    return result;
  };

  struct Details {};
  Details getCurrentLayerDetails() { return Details{}; }
  static RealT strengthHook(Details details,
                            RealT strength,
                            RealT localSlipRate,
                            RealT deltaT,
                            RealT vStar,
                            RealT prakashLength,
                            size_t ltsFace,
                            size_t pointIndex) {
    return strength;
  };
};

template <typename Config>
class BiMaterialFault {
  public:
  using RealT = typename Config::RealT;
  BiMaterialFault(DRParameters* parameters){};

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture const* const dynRup,
                          RealT fullUpdateTime) {
    auto* concreteLts =
        dynamic_cast<seissol::initializers::LTSLinearSlipWeakeningBimaterial const* const>(dynRup);
    this->regularisedStrength = layerData.var(concreteLts->regularisedStrength);
  }

  static RealT resampleSlipRate([[maybe_unused]] RealT const* resampleMatrix,
                                RealT const (&slipRateMagnitude)[dr::misc::numPaddedPoints<Config>],
                                size_t pointIndex) {
    return slipRateMagnitude[pointIndex];
  };

  struct Details {
    RealT (*regularisedStrength)[misc::numPaddedPoints<Config>];
  };

  Details getCurrentLayerDetails() {
    Details details{this->regularisedStrength};
    return details;
  }

  static RealT strengthHook(Details details,
                            RealT faultStrength,
                            RealT localSlipRate,
                            RealT deltaT,
                            RealT vStar,
                            RealT prakashLength,
                            size_t ltsFace,
                            size_t pointIndex) {

    auto* regularisedStrength = details.regularisedStrength[ltsFace];
    assert(regularisedStrength != nullptr && "regularisedStrength is not initialized");

    const RealT expterm = sycl::exp(-(sycl::max(static_cast<RealT>(0.0), localSlipRate) + vStar) *
                                    deltaT / prakashLength);

    const RealT newStrength =
        regularisedStrength[pointIndex] * expterm + faultStrength * (1.0 - expterm);

    regularisedStrength[pointIndex] = newStrength;
    return newStrength;
  };

  private:
  RealT (*regularisedStrength)[misc::numPaddedPoints<Config>];
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPU_LINEARSLIPWEAKENING_H
