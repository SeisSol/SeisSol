#ifndef SEISSOL_GPU_RATEANDSTATEFASTVELOCITYWEAKENING_H
#define SEISSOL_GPU_RATEANDSTATEFASTVELOCITYWEAKENING_H

#include "DynamicRupture/FrictionLaws/GpuImpl/RateAndState.h"

namespace seissol::dr::friction_law::gpu {

template <typename TPMethod>
class FastVelocityWeakeningLaw
    : public RateAndStateBase<FastVelocityWeakeningLaw<TPMethod>, TPMethod> {
  public:
  using RateAndStateBase<FastVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime) {}

  void copySpecificLtsDataTreeToLocal(seissol::initializer::Layer& layerData,
                                      const seissol::initializer::DynamicRupture* const dynRup,
                                      real fullUpdateTime) override {
    using SelfInitializerType = seissol::initializer::LTSRateAndStateFastVelocityWeakening;
    auto* concreteLts = dynamic_cast<const SelfInitializerType*>(dynRup);
    this->srW = layerData.var(concreteLts->rsSrW, seissol::initializer::AllocationPlace::Device);

    using ParentType = RateAndStateBase<FastVelocityWeakeningLaw<TPMethod>, TPMethod>;
    ParentType::copySpecificLtsDataTreeToLocal(layerData, dynRup, fullUpdateTime);
  }

  struct Details {
    decltype(FastVelocityWeakeningLaw::a) a;
    decltype(FastVelocityWeakeningLaw::sl0) sl0;
    decltype(FastVelocityWeakeningLaw::srW) srW;
    decltype(seissol::initializer::parameters::DRParameters::rsSr0) rsSr0;
    decltype(seissol::initializer::parameters::DRParameters::rsF0) rsF0;
    decltype(seissol::initializer::parameters::DRParameters::rsB) rsB;
  };

  Details getCurrentLtsLayerDetails() {
    Details details{};
    details.a = this->a;
    details.sl0 = this->sl0;
    details.srW = this->srW;
    details.rsSr0 = this->drParameters->rsSr0;
    details.rsF0 = this->drParameters->rsF0;
    details.rsB = this->drParameters->rsB;
    return details;
  }

  void updateStateVariable(double timeIncrement) {
    auto* devStateVarReference{this->initialVariables.stateVarReference};
    auto* devLocalSlipRate{this->initialVariables.localSlipRate};
    auto* devStateVariableBuffer{this->stateVariableBuffer};
    const double muW{this->drParameters->muW};
    auto details = this->getCurrentLtsLayerDetails();

    sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints}, {misc::NumPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        const double localSl0 = details.sl0[ltsFace][pointIndex];
        const double localA = details.a[ltsFace][pointIndex];
        const double localSrW = details.srW[ltsFace][pointIndex];
        const double localSlipRate = devLocalSlipRate[ltsFace][pointIndex];

        const double lowVelocityFriction =
            details.rsF0 - (details.rsB - localA) * sycl::log(localSlipRate / details.rsSr0);

        const double steadyStateFrictionCoefficient =
            muW + (lowVelocityFriction - muW) /
                      sycl::pow(1.0 + sycl::pown(localSlipRate / localSrW, 8), 1.0 / 8.0);

        const double steadyStateStateVariable =
            localA * sycl::log(details.rsSr0 / localSlipRate *
                               (sycl::exp(steadyStateFrictionCoefficient / localA) -
                                sycl::exp(-steadyStateFrictionCoefficient / localA)));

        const double exp1 = sycl::exp(-localSlipRate * (timeIncrement / localSl0));
        const double localStateVariable = steadyStateStateVariable * (1.0 - exp1) +
                                          exp1 * devStateVarReference[ltsFace][pointIndex];

        devStateVariableBuffer[ltsFace][pointIndex] = localStateVariable;
      });
    });
  }

  static double updateMu(double localSlipRateMagnitude,
                         double localStateVariable,
                         Details details,
                         size_t ltsFace,
                         size_t pointIndex) {
    const double localA = details.a[ltsFace][pointIndex];
    const double x =
        0.5 / details.rsSr0 * sycl::exp(localStateVariable / localA) * localSlipRateMagnitude;
    return localA * sycl::asinh(x);
  }

  static double updateMuDerivative(double localSlipRateMagnitude,
                                   double localStateVariable,
                                   Details details,
                                   size_t ltsFace,
                                   size_t pointIndex) {
    const double localA = details.a[ltsFace][pointIndex];
    const double c = 0.5 / details.rsSr0 * sycl::exp(localStateVariable / localA);
    return localA * c / sycl::sqrt(sycl::pown(localSlipRateMagnitude * c, 2) + 1.0);
  }

  void resampleStateVar(real (*devStateVariableBuffer)[misc::NumPaddedPoints]) {
    auto* devStateVariable{this->stateVariable};
    auto* resampleMatrix{this->resampleMatrix};

    constexpr auto Dim0 = misc::dimSize<init::resample, 0>();
    constexpr auto Dim1 = misc::dimSize<init::resample, 1>();
    static_assert(Dim0 == misc::NumPaddedPoints);
    static_assert(Dim0 >= Dim1);

    sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints}, {misc::NumPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      // NOLINTNEXTLINE
      sycl::accessor<real, 1, sycl::access::mode::read_write, sycl::access::target::local>
          deltaStateVar(misc::NumPaddedPoints, cgh);

      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        const auto localStateVariable = devStateVariable[ltsFace][pointIndex];
        deltaStateVar[pointIndex] =
            devStateVariableBuffer[ltsFace][pointIndex] - localStateVariable;
        item.barrier(sycl::access::fence_space::local_space);

        real resampledDeltaStateVar{0.0};
        for (size_t i{0}; i < Dim1; ++i) {
          resampledDeltaStateVar += resampleMatrix[pointIndex + i * Dim0] * deltaStateVar[i];
        }

        devStateVariable[ltsFace][pointIndex] = localStateVariable + resampledDeltaStateVar;
      });
    });
  }

  void executeIfNotConverged() {}

  protected:
  real (*srW)[misc::NumPaddedPoints];
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPU_RATEANDSTATEFASTVELOCITYWEAKENING_H
