#ifndef SEISSOL_GPU_RATEANDSTATEFASTVELOCITYWEAKENING_H
#define SEISSOL_GPU_RATEANDSTATEFASTVELOCITYWEAKENING_H

#include "DynamicRupture/FrictionLaws/GpuImpl/RateAndState.h"

namespace seissol::dr::friction_law::gpu {

template <typename TPMethod>
class FastVelocityWeakeningLaw
    : public RateAndStateBase<FastVelocityWeakeningLaw<TPMethod>, TPMethod> {
  public:
  using RateAndStateBase<FastVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture const* const dynRup,
                          real fullUpdateTime) {}

  void copySpecificLtsDataTreeToLocal(seissol::initializers::Layer& layerData,
                                      seissol::initializers::DynamicRupture const* const dynRup,
                                      real fullUpdateTime) override {
    using SelfInitializerType = seissol::initializers::LTSRateAndStateFastVelocityWeakening;
    auto* concreteLts = dynamic_cast<SelfInitializerType const* const>(dynRup);
    this->srW = layerData.var(concreteLts->rsSrW);

    using ParentType = RateAndStateBase<FastVelocityWeakeningLaw<TPMethod>, TPMethod>;
    ParentType::copySpecificLtsDataTreeToLocal(layerData, dynRup, fullUpdateTime);
  }

  #pragma omp declare target
  struct Details {
    decltype(FastVelocityWeakeningLaw::a) a;
    decltype(FastVelocityWeakeningLaw::sl0) sl0;
    decltype(FastVelocityWeakeningLaw::srW) srW;
    decltype(dr::DRParameters::rsSr0) rsSr0;
    decltype(dr::DRParameters::rsF0) rsF0;
    decltype(dr::DRParameters::rsB) rsB;
    std::size_t layerSize;
  };

  #pragma omp declare mapper(Details det) map(det, det.a[0:det.layerSize], det.srW[0:det.layerSize], det.sl0[0:det.layerSize])

  Details getCurrentLtsLayerDetails() {
    Details details{};
    details.a = this->a;
    details.sl0 = this->sl0;
    details.srW = this->srW;
    details.rsSr0 = this->drParameters->rsSr0;
    details.rsF0 = this->drParameters->rsF0;
    details.rsB = this->drParameters->rsB;
    details.layerSize = this->currLayerSize;
    return details;
  }

  void updateStateVariable(double timeIncrement) {
    const auto layerSize{this->currLayerSize};
    auto* devStateVarReference{this->initialVariables.stateVarReference};
    auto* devLocalSlipRate{this->initialVariables.localSlipRate};
    auto* devStateVariableBuffer{this->stateVariableBuffer};
    const double muW{this->drParameters->muW};
    auto details = this->getCurrentLtsLayerDetails();

    auto* queue{this->queue};

    // #pragma omp distribute
    #pragma omp target teams distribute depend(inout: *queue) device(TARGETDART_ANY) map(to: details, details.sl0[0:layerSize], details.a[0:layerSize], details.srW[0:layerSize], devStateVarReference[0:layerSize], devLocalSlipRate[0:layerSize]) map(from: devStateVariableBuffer[0:layerSize]) nowait
      for (int ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
        #pragma omp parallel for schedule(static, 1)
        for (int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {

        const double localSl0 = details.sl0[ltsFace][pointIndex];
        const double localA = details.a[ltsFace][pointIndex];
        const double localSrW = details.srW[ltsFace][pointIndex];
        const double localSlipRate = devLocalSlipRate[ltsFace][pointIndex];

        const double lowVelocityFriction =
            details.rsF0 - (details.rsB - localA) * std::log(localSlipRate / details.rsSr0);

        const double steadyStateFrictionCoefficient =
            muW + (lowVelocityFriction - muW) /
                      std::pow(1.0 + std::pow(localSlipRate / localSrW, 8), 1.0 / 8.0);

        const double steadyStateStateVariable =
            localA * std::log(details.rsSr0 / localSlipRate *
                               (std::exp(steadyStateFrictionCoefficient / localA) -
                                std::exp(-steadyStateFrictionCoefficient / localA)));

        const double exp1 = std::exp(-localSlipRate * (timeIncrement / localSl0));
        const double localStateVariable = steadyStateStateVariable * (1.0 - exp1) +
                                          exp1 * devStateVarReference[ltsFace][pointIndex];

        devStateVariableBuffer[ltsFace][pointIndex] = localStateVariable;
      }
    }
  }

  static double updateMu(double localSlipRateMagnitude,
                         double localStateVariable,
                         Details details,
                         size_t ltsFace,
                         size_t pointIndex) {
    const double localA = details.a[ltsFace][pointIndex];
    const double x =
        0.5 / details.rsSr0 * std::exp(localStateVariable / localA) * localSlipRateMagnitude;
    return localA * std::asinh(x);
  }

  static double updateMuDerivative(double localSlipRateMagnitude,
                                   double localStateVariable,
                                   Details details,
                                   size_t ltsFace,
                                   size_t pointIndex) {
    const double localA = details.a[ltsFace][pointIndex];
    const double c = 0.5 / details.rsSr0 * std::exp(localStateVariable / localA);
    return localA * c / std::sqrt(std::pow(localSlipRateMagnitude * c, 2) + 1.0);
  }

  void resampleStateVar(real (*devStateVariableBuffer)[misc::numPaddedPoints]) {
    auto* devStateVariable{this->stateVariable};
    auto* resampleMatrix{this->resampleMatrix};

    const auto layerSize{this->currLayerSize};
    constexpr auto dim0 = misc::dimSize<init::resample, 0>();
    constexpr auto dim1 = misc::dimSize<init::resample, 1>();
    constexpr auto resampleSize = dim0 * dim1 * sizeof(real);
    static_assert(dim0 == misc::numPaddedPoints);
    static_assert(dim0 >= dim1);

    auto* queue{this->queue};

     /* std::accessor<real, 1, std::access::mode::read_write, std::access::target::local>
          deltaStateVar(misc::numPaddedPoints, cgh);*/
    // #pragma omp distribute
    #pragma omp target teams distribute depend(inout: *queue) device(TARGETDART_ANY) map(to: devStateVariableBuffer[0:layerSize], resampleMatrix[0:resampleSize]) map(tofrom: devStateVariable[0:layerSize]) nowait
    // allocate(omp_pteam_mem_alloc:deltaStateVar)
      for (int ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
        real deltaStateVar[misc::numPaddedPoints];
        #pragma omp parallel for schedule(static, 1)
        for (int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
          deltaStateVar[pointIndex] =
              devStateVariableBuffer[ltsFace][pointIndex] - devStateVariable[ltsFace][pointIndex];
        }
        #pragma omp parallel for schedule(static, 1)
        for (int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        real resampledDeltaStateVar{0.0};
        for (size_t i{0}; i < dim1; ++i) {
          resampledDeltaStateVar += resampleMatrix[pointIndex + i * dim0] * deltaStateVar[i];
        }

        devStateVariable[ltsFace][pointIndex] += resampledDeltaStateVar;
      }
    }
  }

  void executeIfNotConverged() {}

  protected:
  real (*srW)[misc::numPaddedPoints];
  #pragma omp end declare target
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPU_RATEANDSTATEFASTVELOCITYWEAKENING_H
