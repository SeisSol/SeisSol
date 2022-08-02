#include "Factory.h"

#include "FrictionLaws/FrictionLaws.h"
#include "FrictionLaws/ThermalPressurization/NoTP.h"
#include "FrictionLaws/ThermalPressurization/ThermalPressurization.h"

#ifdef ACL_DEVICE_OFFLOAD
namespace friction_law_impl = seissol::dr::friction_law::gpu;
#else
namespace friction_law_impl = seissol::dr::friction_law;
#endif

namespace seissol::dr::factory {
std::unique_ptr<AbstractFactory> getFactory(std::shared_ptr<dr::DRParameters> drParameters) {
  switch (drParameters->frictionLawType) {
  case FrictionLawType::NoFault:
    return std::make_unique<NoFaultFactory>(drParameters);
  case FrictionLawType::ImposedSlipRatesYoffe:
    return std::make_unique<ImposedSlipRatesYoffeFactory>(drParameters);
  case FrictionLawType::ImposedSlipRatesGaussian:
    return std::make_unique<ImposedSlipRatesGaussianFactory>(drParameters);
  case FrictionLawType::LinearSlipWeakening:
    return std::make_unique<LinearSlipWeakeningFactory>(drParameters);
  case FrictionLawType::LinearSlipWeakeningBimaterial:
    return std::make_unique<LinearSlipWeakeningBimaterialFactory>(drParameters);
  case FrictionLawType::RateAndStateAgingLaw:
    return std::make_unique<RateAndStateAgingFactory>(drParameters);
  case FrictionLawType::RateAndStateSlipLaw:
    return std::make_unique<RateAndStateSlipFactory>(drParameters);
  case FrictionLawType::RateAndStateVelocityWeakening:
    logError() << "friction law 7 currently disabled";
    return std::unique_ptr<AbstractFactory>(nullptr);
  case FrictionLawType::RateAndStateAgingNucleation:
    logError() << "friction law 101 currently disabled";
    return std::unique_ptr<AbstractFactory>(nullptr);
  case FrictionLawType::RateAndStateFastVelocityWeakening:
    return std::make_unique<RateAndStateFastVelocityWeakeningFactory>(drParameters);
  default:
    logError() << "unknown friction law";
    return nullptr;
  }
}

DynamicRuptureTuple NoFaultFactory::produce() {
  return {std::make_unique<seissol::initializers::DynamicRupture>(),
          std::make_unique<initializers::NoFaultInitializer>(drParameters),
          std::make_unique<friction_law::NoFault>(drParameters.get()),
          std::make_unique<output::OutputManager>(new output::NoFault)};
}

DynamicRuptureTuple LinearSlipWeakeningFactory::produce() {
  return {std::make_unique<seissol::initializers::LTS_LinearSlipWeakening>(),
          std::make_unique<initializers::LinearSlipWeakeningInitializer>(drParameters),
          std::make_unique<
              friction_law_impl::LinearSlipWeakeningLaw<friction_law_impl::NoSpecialization>>(
              drParameters.get()),
          std::make_unique<output::OutputManager>(new output::LinearSlipWeakening)};
}

DynamicRuptureTuple RateAndStateAgingFactory::produce() {
  if (drParameters->isThermalPressureOn) {
    return {std::make_unique<seissol::initializers::LTS_RateAndState>(),
            std::make_unique<initializers::RateAndStateInitializer>(drParameters),
            std::make_unique<friction_law::AgingLaw<friction_law::ThermalPressurization>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(new output::RateAndStateThermalPressurization)};
  } else {
    return {std::make_unique<seissol::initializers::LTS_RateAndState>(),
            std::make_unique<initializers::RateAndStateInitializer>(drParameters),
            std::make_unique<friction_law::AgingLaw<friction_law::NoTP>>(drParameters.get()),
            std::make_unique<output::OutputManager>(new output::RateAndState)};
  }
}

DynamicRuptureTuple RateAndStateSlipFactory::produce() {
  if (drParameters->isThermalPressureOn) {
    return {std::make_unique<seissol::initializers::LTS_RateAndState>(),
            std::make_unique<initializers::RateAndStateInitializer>(drParameters),
            std::make_unique<friction_law::SlipLaw<friction_law::ThermalPressurization>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(new output::RateAndStateThermalPressurization)};
  } else {
    return {std::make_unique<seissol::initializers::LTS_RateAndState>(),
            std::make_unique<initializers::RateAndStateInitializer>(drParameters),
            std::make_unique<friction_law::SlipLaw<friction_law::NoTP>>(drParameters.get()),
            std::make_unique<output::OutputManager>(new output::RateAndState)};
  }
}

DynamicRuptureTuple LinearSlipWeakeningBimaterialFactory::produce() {
  return {std::make_unique<seissol::initializers::LTS_LinearSlipWeakeningBimaterial>(),
          std::make_unique<initializers::LinearSlipWeakeningBimaterialInitializer>(drParameters),
          std::make_unique<friction_law::LinearSlipWeakeningLaw<friction_law::BiMaterialFault>>(
              drParameters.get()),
          std::make_unique<output::OutputManager>(new output::LinearSlipWeakeningBimaterial)};
}

DynamicRuptureTuple ImposedSlipRatesYoffeFactory::produce() {
  return {
      std::make_unique<seissol::initializers::LTS_ImposedSlipRatesYoffe>(),
      std::make_unique<initializers::ImposedSlipRatesYoffeInitializer>(drParameters),
      std::make_unique<friction_law::ImposedSlipRates<friction_law::YoffeSTF>>(drParameters.get()),
      std::make_unique<output::OutputManager>(new output::ImposedSlipRates)};
}

DynamicRuptureTuple ImposedSlipRatesGaussianFactory::produce() {
  return {std::make_unique<seissol::initializers::LTS_ImposedSlipRatesGaussian>(),
          std::make_unique<initializers::ImposedSlipRatesGaussianInitializer>(drParameters),
          std::make_unique<friction_law::ImposedSlipRates<friction_law::GaussianSTF>>(
              drParameters.get()),
          std::make_unique<output::OutputManager>(new output::ImposedSlipRates)};
}

DynamicRuptureTuple RateAndStateFastVelocityWeakeningFactory::produce() {
  if (drParameters->isThermalPressureOn) {
    return {
        std::make_unique<seissol::initializers::LTS_RateAndStateThermalPressurization>(),
        std::make_unique<initializers::RateAndStateThermalPressurizationInitializer>(drParameters),
        std::make_unique<
            friction_law::FastVelocityWeakeningLaw<friction_law::ThermalPressurization>>(
            drParameters.get()),
        std::make_unique<output::OutputManager>(new output::RateAndStateThermalPressurization)};
  } else {
    return {std::make_unique<seissol::initializers::LTS_RateAndStateFastVelocityWeakening>(),
            std::make_unique<initializers::RateAndStateFastVelocityInitializer>(drParameters),
            std::make_unique<friction_law::FastVelocityWeakeningLaw<friction_law::NoTP>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(new output::RateAndState)};
  }
}
} // namespace seissol::dr::factory
