#include "Factory.h"

#include "FrictionLaws/FrictionLaws.h"
#include "FrictionLaws/ThermalPressurization/ThermalPressurization.h"

#ifdef ACL_DEVICE
namespace friction_law_impl = seissol::dr::friction_law::gpu;
using MathFunctions = seissol::functions::SyclStdFunctions;
#else
namespace friction_law_impl = seissol::dr::friction_law;
using MathFunctions = seissol::functions::HostStdFunctions;
#endif

namespace seissol::dr::factory {
std::unique_ptr<AbstractFactory>
    getFactory(std::shared_ptr<seissol::initializer::parameters::DRParameters> drParameters,
               seissol::SeisSol& seissolInstance) {
  switch (drParameters->frictionLawType) {
  case seissol::initializer::parameters::FrictionLawType::NoFault:
    return std::make_unique<NoFaultFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::ImposedSlipRatesYoffe:
    return std::make_unique<ImposedSlipRatesYoffeFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::ImposedSlipRatesGaussian:
    return std::make_unique<ImposedSlipRatesGaussianFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::LinearSlipWeakening:
    return std::make_unique<LinearSlipWeakeningFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::LinearSlipWeakeningBimaterial:
    return std::make_unique<LinearSlipWeakeningBimaterialFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateAgingLaw:
    return std::make_unique<RateAndStateAgingFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateSlipLaw:
    return std::make_unique<RateAndStateSlipFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateVelocityWeakening:
    logError() << "friction law 7 currently disabled";
    return std::unique_ptr<AbstractFactory>(nullptr);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateAgingNucleation:
    logError() << "friction law 101 currently disabled";
    return std::unique_ptr<AbstractFactory>(nullptr);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateFastVelocityWeakening:
    return std::make_unique<RateAndStateFastVelocityWeakeningFactory>(drParameters,
                                                                      seissolInstance);
  default:
    logError() << "unknown friction law";
    return nullptr;
  }
}

DynamicRuptureTuple NoFaultFactory::produce() {
  return {std::make_unique<seissol::initializer::DynamicRupture>(),
          std::make_unique<initializer::NoFaultInitializer>(drParameters, seissolInstance),
          std::make_unique<friction_law::NoFault>(drParameters.get()),
          std::make_unique<output::OutputManager>(std::make_unique<output::NoFault>(),
                                                  seissolInstance)};
}

DynamicRuptureTuple LinearSlipWeakeningFactory::produce() {
  return {
      std::make_unique<seissol::initializer::LTSLinearSlipWeakening>(),
      std::make_unique<initializer::LinearSlipWeakeningInitializer>(drParameters, seissolInstance),
      std::make_unique<
          friction_law_impl::LinearSlipWeakeningLaw<friction_law_impl::NoSpecialization>>(
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
            std::make_unique<output::OutputManager>(
                std::make_unique<output::RateAndStateThermalPressurization>(), seissolInstance)};
  } else {
    return {
        std::make_unique<seissol::initializer::LTSRateAndState>(),
        std::make_unique<initializer::RateAndStateInitializer>(drParameters, seissolInstance),
        std::make_unique<friction_law_impl::AgingLaw<friction_law_impl::NoTP>>(drParameters.get()),
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
            std::make_unique<output::OutputManager>(
                std::make_unique<output::RateAndStateThermalPressurization>(), seissolInstance)};
  } else {
    return {
        std::make_unique<seissol::initializer::LTSRateAndState>(),
        std::make_unique<initializer::RateAndStateInitializer>(drParameters, seissolInstance),
        std::make_unique<friction_law_impl::SlipLaw<friction_law_impl::NoTP>>(drParameters.get()),
        std::make_unique<output::OutputManager>(std::make_unique<output::RateAndState>(),
                                                seissolInstance)};
  }
}

DynamicRuptureTuple LinearSlipWeakeningBimaterialFactory::produce() {
  using Specialization = friction_law_impl::BiMaterialFault;
  using FrictionLawType = friction_law_impl::LinearSlipWeakeningLaw<Specialization>;

  return {std::make_unique<seissol::initializer::LTSLinearSlipWeakeningBimaterial>(),
          std::make_unique<initializer::LinearSlipWeakeningBimaterialInitializer>(drParameters,
                                                                                  seissolInstance),
          std::make_unique<FrictionLawType>(drParameters.get()),
          std::make_unique<output::OutputManager>(
              std::make_unique<output::LinearSlipWeakeningBimaterial>(), seissolInstance)};
}

DynamicRuptureTuple ImposedSlipRatesYoffeFactory::produce() {
  return {
      std::make_unique<seissol::initializer::LTSImposedSlipRatesYoffe>(),
      std::make_unique<initializer::ImposedSlipRatesYoffeInitializer>(drParameters,
                                                                      seissolInstance),
      std::make_unique<friction_law_impl::ImposedSlipRates<friction_law::YoffeSTF<MathFunctions>>>(
          drParameters.get()),
      std::make_unique<output::OutputManager>(std::make_unique<output::ImposedSlipRates>(),
                                              seissolInstance)};
}

DynamicRuptureTuple ImposedSlipRatesGaussianFactory::produce() {
  return {std::make_unique<seissol::initializer::LTSImposedSlipRatesGaussian>(),
          std::make_unique<initializer::ImposedSlipRatesGaussianInitializer>(drParameters,
                                                                             seissolInstance),
          std::make_unique<
              friction_law_impl::ImposedSlipRates<friction_law::GaussianSTF<MathFunctions>>>(
              drParameters.get()),
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
            std::make_unique<output::OutputManager>(
                std::make_unique<output::RateAndStateThermalPressurization>(), seissolInstance)};
  } else {
    return {std::make_unique<seissol::initializer::LTSRateAndStateFastVelocityWeakening>(),
            std::make_unique<initializer::RateAndStateFastVelocityInitializer>(drParameters,
                                                                               seissolInstance),
            std::make_unique<friction_law_impl::FastVelocityWeakeningLaw<friction_law_impl::NoTP>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(std::make_unique<output::RateAndState>(),
                                                    seissolInstance)};
  }
}
} // namespace seissol::dr::factory
