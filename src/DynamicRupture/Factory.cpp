#include "Factory.h"

#include "FrictionLaws/FrictionLaws.h"
#include "FrictionLaws/ThermalPressurization/NoTP.h"
#include "FrictionLaws/ThermalPressurization/ThermalPressurization.h"
#include <Solver/Interoperability.h>

namespace seissol::dr::factory {
std::unique_ptr<AbstractFactory> getFactory(dr::DRParameters& drParameters) {
  switch (drParameters.frictionLawType) {
  case FrictionLawType::NoFault:
    return std::make_unique<NoFaultFactory>(drParameters);
  case FrictionLawType::ImposedSlipRates:
    return std::make_unique<ImposedSlipRatesFactory>(drParameters);
  case FrictionLawType::LinearSlipWeakening:
    return std::make_unique<LinearSlipWeakeningFactory>(drParameters);
  case FrictionLawType::LinearSlipWeakeningForcedRuptureTime:
    return std::make_unique<LinearSlipWeakeningForcedRuptureTimeFactory>(drParameters);
  // Prakash-Clifton regularisation for bimaterial faults: see (Pelties et al. 2014)
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

Products NoFaultFactory::produce() {
  return {std::make_unique<seissol::initializers::DynamicRupture>(),
          std::make_unique<initializers::NoFaultInitializer>(drParameters),
          std::make_unique<friction_law::NoFault>(drParameters),
          std::make_unique<output::NoFault>()};
}

Products LinearSlipWeakeningFactory::produce() {
  return {std::make_unique<seissol::initializers::LTS_LinearSlipWeakening>(),
          std::make_unique<initializers::LinearSlipWeakeningInitializer>(drParameters),
          std::make_unique<friction_law::LinearSlipWeakeningLaw>(drParameters),
          std::make_unique<output::LinearSlipWeakening>()};
}

Products RateAndStateAgingFactory::produce() {
  if (drParameters.isThermalPressureOn) {
    return {
        std::make_unique<seissol::initializers::LTS_RateAndState>(),
        std::make_unique<initializers::RateAndStateInitializer>(drParameters),
        std::make_unique<friction_law::AgingLaw<friction_law::ThermalPressurization>>(drParameters),
        std::make_unique<output::RateAndState>()};
  } else {
    return {std::make_unique<seissol::initializers::LTS_RateAndState>(),
            std::make_unique<initializers::RateAndStateInitializer>(drParameters),
            std::make_unique<friction_law::AgingLaw<friction_law::NoTP>>(drParameters),
            std::make_unique<output::RateAndState>()};
  }
}

Products RateAndStateSlipFactory::produce() {
  if (drParameters.isThermalPressureOn) {
    return {
        std::make_unique<seissol::initializers::LTS_RateAndState>(),
        std::make_unique<initializers::RateAndStateInitializer>(drParameters),
        std::make_unique<friction_law::SlipLaw<friction_law::ThermalPressurization>>(drParameters),
        std::make_unique<output::RateAndState>()};
  } else {
    return {std::make_unique<seissol::initializers::LTS_RateAndState>(),
            std::make_unique<initializers::RateAndStateInitializer>(drParameters),
            std::make_unique<friction_law::SlipLaw<friction_law::NoTP>>(drParameters),
            std::make_unique<output::RateAndState>()};
  }
}

Products LinearSlipWeakeningBimaterialFactory::produce() {
  return {std::make_unique<seissol::initializers::LTS_LinearSlipWeakeningBimaterial>(),
          std::make_unique<initializers::LinearSlipWeakeningBimaterialInitializer>(drParameters),
          std::make_unique<friction_law::LinearSlipWeakeningLawBimaterial>(drParameters),
          std::make_unique<output::LinearSlipWeakeningBimaterial>()};
}

Products LinearSlipWeakeningForcedRuptureTimeFactory::produce() {
  return {
      std::make_unique<seissol::initializers::LTS_LinearSlipWeakeningForcedRuptureTime>(),
      std::make_unique<initializers::LinearSlipWeakeningForcedRuptureTimeInitializer>(drParameters),
      std::make_unique<friction_law::LinearSlipWeakeningLawForcedRuptureTime>(drParameters),
      std::make_unique<output::LinearSlipWeakening>()};
}

Products ImposedSlipRatesFactory::produce() {
  return {std::make_unique<seissol::initializers::LTS_ImposedSlipRates>(),
          std::make_unique<initializers::ImposedSlipRatesInitializer>(drParameters),
          std::make_unique<friction_law::ImposedSlipRates>(drParameters),
          std::make_unique<output::ImposedSlipRates>()};
}

Products RateAndStateFastVelocityWeakeningFactory::produce() {
  if (drParameters.isThermalPressureOn) {
    return {std::make_unique<seissol::initializers::LTS_RateAndState>(),
            std::make_unique<initializers::RateAndStateInitializer>(drParameters),
            std::make_unique<
                friction_law::FastVelocityWeakeningLaw<friction_law::ThermalPressurization>>(
                drParameters),
            std::make_unique<output::RateAndState>()};
  } else {
    return {
        std::make_unique<seissol::initializers::LTS_RateAndState>(),
        std::make_unique<initializers::RateAndStateInitializer>(drParameters),
        std::make_unique<friction_law::FastVelocityWeakeningLaw<friction_law::NoTP>>(drParameters),
        std::make_unique<output::RateAndState>()};
  }
}
} // namespace seissol::dr::factory
