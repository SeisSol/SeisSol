#include "Factory.h"

#include "FrictionLaws/FrictionLaws.h"
#include "FrictionLaws/ThermalPressurization/ThermalPressurization.h"
#include "Initializer/DynamicRupture.h"
#include "Initializer/Initializers.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Output/Output.h"
#include <memory>
#include <utils/logger.h>

// for now, fake the friction laws to be on CPU instead, if we have no GPU
#ifdef ACL_DEVICE
namespace friction_law_gpu = seissol::dr::friction_law::gpu;
#else
namespace friction_law_gpu = seissol::dr::friction_law;
#endif

namespace seissol::dr::factory {
std::unique_ptr<AbstractFactory>
    getFactory(const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
               seissol::SeisSol& seissolInstance) {
  switch (drParameters->frictionLawType) {
  case seissol::initializer::parameters::FrictionLawType::NoFault:
    return std::make_unique<NoFaultFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::ImposedSlipRatesYoffe:
    return std::make_unique<ImposedSlipRatesYoffeFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::ImposedSlipRatesGaussian:
    return std::make_unique<ImposedSlipRatesGaussianFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::ImposedSlipRatesDelta:
    return std::make_unique<ImposedSlipRatesDeltaFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::LinearSlipWeakening:
    return std::make_unique<LinearSlipWeakeningFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::LinearSlipWeakeningBimaterial:
    return std::make_unique<LinearSlipWeakeningBimaterialFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::LinearSlipWeakeningTPApprox:
    return std::make_unique<LinearSlipWeakeningTPApproxFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateAgingLaw:
    return std::make_unique<RateAndStateAgingFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateSlipLaw:
    return std::make_unique<RateAndStateSlipFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateVelocityWeakening:
    logError() << "friction law 7 currently disabled";
    return {nullptr};
  case seissol::initializer::parameters::FrictionLawType::RateAndStateAgingNucleation:
    logError() << "friction law 101 currently disabled";
    return {nullptr};
  case seissol::initializer::parameters::FrictionLawType::RateAndStateFastVelocityWeakening:
    return std::make_unique<RateAndStateFastVelocityWeakeningFactory>(drParameters,
                                                                      seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::AdjointSlip:
    return std::make_unique<AdjointRSFSlipFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::AdjointRSFFastVelWeakening:
    return std::make_unique<AdjointRSFFastVelWeakeningFactory>(drParameters, seissolInstance);
  default:
    logError() << "unknown friction law";
    return nullptr;
  }
}

DynamicRuptureTuple NoFaultFactory::produce() {
  return {std::make_unique<seissol::initializer::DynamicRupture>(),
          std::make_unique<initializer::NoFaultInitializer>(drParameters, seissolInstance),
          std::make_unique<friction_law::NoFault>(drParameters.get()),
          std::make_unique<friction_law_gpu::NoFault>(drParameters.get()),
          std::make_unique<output::OutputManager>(std::make_unique<output::NoFault>(),
                                                  seissolInstance)};
}

DynamicRuptureTuple LinearSlipWeakeningFactory::produce() {
  return {
      std::make_unique<seissol::initializer::LTSLinearSlipWeakening>(),
      std::make_unique<initializer::LinearSlipWeakeningInitializer>(drParameters, seissolInstance),
      std::make_unique<friction_law::LinearSlipWeakeningLaw<friction_law::NoSpecialization>>(
          drParameters.get()),
      std::make_unique<
          friction_law_gpu::LinearSlipWeakeningLaw<friction_law_gpu::NoSpecialization>>(
          drParameters.get()),
      std::make_unique<output::OutputManager>(std::make_unique<output::LinearSlipWeakening>(),
                                              seissolInstance)};
}

DynamicRuptureTuple RateAndStateAgingFactory::produce() {
  if (drParameters->isThermalPressureOn) {
    return {std::make_unique<seissol::initializer::LTSRateAndState>(),
            std::make_unique<initializer::RateAndStateInitializer>(drParameters, seissolInstance),
            std::make_unique<friction_law::AgingLaw<friction_law::ThermalPressurization>>(
                drParameters.get()),
            std::make_unique<friction_law::AgingLaw<friction_law::ThermalPressurization>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(
                std::make_unique<output::RateAndStateThermalPressurization>(), seissolInstance)};
  } else {
    return {
        std::make_unique<seissol::initializer::LTSRateAndState>(),
        std::make_unique<initializer::RateAndStateInitializer>(drParameters, seissolInstance),
        std::make_unique<friction_law::AgingLaw<friction_law::NoTP>>(drParameters.get()),
        std::make_unique<friction_law_gpu::AgingLaw<friction_law_gpu::NoTP>>(drParameters.get()),
        std::make_unique<output::OutputManager>(std::make_unique<output::RateAndState>(),
                                                seissolInstance)};
  }
}

DynamicRuptureTuple RateAndStateSlipFactory::produce() {
  if (drParameters->isThermalPressureOn) {
    return {std::make_unique<seissol::initializer::LTSRateAndState>(),
            std::make_unique<initializer::RateAndStateInitializer>(drParameters, seissolInstance),
            std::make_unique<friction_law::SlipLaw<friction_law::ThermalPressurization>>(
                drParameters.get()),
            std::make_unique<friction_law::SlipLaw<friction_law::ThermalPressurization>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(
                std::make_unique<output::RateAndStateThermalPressurization>(), seissolInstance)};
  } else {
    return {std::make_unique<seissol::initializer::LTSRateAndState>(),
            std::make_unique<initializer::RateAndStateInitializer>(drParameters, seissolInstance),
            std::make_unique<friction_law::SlipLaw<friction_law::NoTP>>(drParameters.get()),
            std::make_unique<friction_law_gpu::SlipLaw<friction_law_gpu::NoTP>>(drParameters.get()),
            std::make_unique<output::OutputManager>(std::make_unique<output::RateAndState>(),
                                                    seissolInstance)};
  }
}

DynamicRuptureTuple LinearSlipWeakeningBimaterialFactory::produce() {
  using Specialization = friction_law::BiMaterialFault;
  using FrictionLawType = friction_law::LinearSlipWeakeningLaw<Specialization>;
  using SpecializationGpu = friction_law_gpu::BiMaterialFault;
  using FrictionLawTypeGpu = friction_law_gpu::LinearSlipWeakeningLaw<SpecializationGpu>;

  return {std::make_unique<seissol::initializer::LTSLinearSlipWeakeningBimaterial>(),
          std::make_unique<initializer::LinearSlipWeakeningBimaterialInitializer>(drParameters,
                                                                                  seissolInstance),
          std::make_unique<FrictionLawType>(drParameters.get()),
          std::make_unique<FrictionLawTypeGpu>(drParameters.get()),
          std::make_unique<output::OutputManager>(
              std::make_unique<output::LinearSlipWeakeningBimaterial>(), seissolInstance)};
}

DynamicRuptureTuple LinearSlipWeakeningTPApproxFactory::produce() {
  using Specialization = friction_law::TPApprox;
  using FrictionLawType = friction_law::LinearSlipWeakeningLaw<Specialization>;
  using SpecializationGpu = friction_law_gpu::TPApprox;
  using FrictionLawTypeGpu = friction_law_gpu::LinearSlipWeakeningLaw<SpecializationGpu>;

  return {
      std::make_unique<seissol::initializer::LTSLinearSlipWeakening>(),
      std::make_unique<initializer::LinearSlipWeakeningInitializer>(drParameters, seissolInstance),
      std::make_unique<FrictionLawType>(drParameters.get()),
      std::make_unique<FrictionLawTypeGpu>(drParameters.get()),
      std::make_unique<output::OutputManager>(std::make_unique<output::LinearSlipWeakening>(),
                                              seissolInstance)};
}

DynamicRuptureTuple ImposedSlipRatesYoffeFactory::produce() {
  return {
      std::make_unique<seissol::initializer::LTSImposedSlipRatesYoffe>(),
      std::make_unique<initializer::ImposedSlipRatesYoffeInitializer>(drParameters,
                                                                      seissolInstance),
      std::make_unique<friction_law::ImposedSlipRates<friction_law::YoffeSTF>>(drParameters.get()),
      std::make_unique<friction_law::ImposedSlipRates<friction_law::YoffeSTF>>(drParameters.get()),
      std::make_unique<output::OutputManager>(std::make_unique<output::ImposedSlipRates>(),
                                              seissolInstance)};
}

DynamicRuptureTuple ImposedSlipRatesGaussianFactory::produce() {
  return {std::make_unique<seissol::initializer::LTSImposedSlipRatesGaussian>(),
          std::make_unique<initializer::ImposedSlipRatesGaussianInitializer>(drParameters,
                                                                             seissolInstance),
          std::make_unique<friction_law::ImposedSlipRates<friction_law::GaussianSTF>>(
              drParameters.get()),
          std::make_unique<friction_law::ImposedSlipRates<friction_law::GaussianSTF>>(
              drParameters.get()),
          std::make_unique<output::OutputManager>(std::make_unique<output::ImposedSlipRates>(),
                                                  seissolInstance)};
}

DynamicRuptureTuple ImposedSlipRatesDeltaFactory::produce() {
  return {
      std::make_unique<seissol::initializer::LTSImposedSlipRatesDelta>(),
      std::make_unique<initializer::ImposedSlipRatesDeltaInitializer>(drParameters,
                                                                      seissolInstance),
      std::make_unique<friction_law::ImposedSlipRates<friction_law::DeltaSTF>>(drParameters.get()),
      std::make_unique<friction_law::ImposedSlipRates<friction_law::DeltaSTF>>(drParameters.get()),
      std::make_unique<output::OutputManager>(std::make_unique<output::ImposedSlipRates>(),
                                              seissolInstance)};
}

DynamicRuptureTuple RateAndStateFastVelocityWeakeningFactory::produce() {
  if (drParameters->isThermalPressureOn) {
    return {std::make_unique<seissol::initializer::LTSRateAndStateThermalPressurization>(),
            std::make_unique<initializer::RateAndStateThermalPressurizationInitializer>(
                drParameters, seissolInstance),
            std::make_unique<
                friction_law::FastVelocityWeakeningLaw<friction_law::ThermalPressurization>>(
                drParameters.get()),
            std::make_unique<
                friction_law::FastVelocityWeakeningLaw<friction_law::ThermalPressurization>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(
                std::make_unique<output::RateAndStateThermalPressurization>(), seissolInstance)};
  } else {
    return {std::make_unique<seissol::initializer::LTSRateAndStateFastVelocityWeakening>(),
            std::make_unique<initializer::RateAndStateFastVelocityInitializer>(drParameters,
                                                                               seissolInstance),
            std::make_unique<friction_law::FastVelocityWeakeningLaw<friction_law::NoTP>>(
                drParameters.get()),
            std::make_unique<friction_law_gpu::FastVelocityWeakeningLaw<friction_law_gpu::NoTP>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(std::make_unique<output::RateAndState>(),
                                                    seissolInstance)};
  }
}

DynamicRuptureTuple AdjointRSFSlipFactory::produce() {
  if (drParameters->isThermalPressureOn) {
    return {std::make_unique<seissol::initializer::LTSAdjointRSF>(),
            std::make_unique<initializer::AdjointRSFInitializer>(drParameters, seissolInstance),
            std::make_unique<friction_law::AdjointSlip<friction_law::ThermalPressurization>>(
                drParameters.get()),
            std::make_unique<friction_law::AdjointSlip<friction_law::ThermalPressurization>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(
                std::make_unique<output::RateAndStateThermalPressurization>(), seissolInstance)};
  } else {
    return {
        std::make_unique<seissol::initializer::LTSAdjointRSF>(),
        std::make_unique<initializer::AdjointRSFInitializer>(drParameters, seissolInstance),
        std::make_unique<friction_law::AdjointSlip<friction_law::NoTP>>(drParameters.get()),
        std::make_unique<friction_law_gpu::AdjointSlip<friction_law_gpu::NoTP>>(drParameters.get()),
        std::make_unique<output::OutputManager>(std::make_unique<output::RateAndState>(),
                                                seissolInstance)};
  }
}

DynamicRuptureTuple AdjointRSFFastVelWeakeningFactory::produce() {
  if (drParameters->isThermalPressureOn) {
    return {std::make_unique<seissol::initializer::LTSRateAndStateThermalPressurization>(),
            std::make_unique<initializer::RateAndStateThermalPressurizationInitializer>(
                drParameters, seissolInstance),
            std::make_unique<
                friction_law::FastVelocityWeakeningLaw<friction_law::ThermalPressurization>>(
                drParameters.get()),
            std::make_unique<
                friction_law::FastVelocityWeakeningLaw<friction_law::ThermalPressurization>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(
                std::make_unique<output::RateAndStateThermalPressurization>(), seissolInstance)};
  } else {
    return {
        std::make_unique<seissol::initializer::LTSAdjointRSFFastVelWeakening>(),
        std::make_unique<initializer::AdjointRSFFastVelInitializer>(drParameters, seissolInstance),
        std::make_unique<friction_law::AdjointFastVelWeakening<friction_law::NoTP>>(
            drParameters.get()),
        std::make_unique<friction_law_gpu::AdjointFastVelWeakening<friction_law_gpu::NoTP>>(
            drParameters.get()),
        std::make_unique<output::OutputManager>(std::make_unique<output::RateAndState>(),
                                                seissolInstance)};
  }
}

} // namespace seissol::dr::factory
