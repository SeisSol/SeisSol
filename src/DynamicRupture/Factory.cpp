#include "Factory.h"

#include <Solver/Interoperability.h>

namespace seissol::dr::factory {
std::unique_ptr<AbstractFactory> getFactory(dr::DRParameters& drParameters) {
  switch (drParameters.frictionLawType) {
  case FrictionLawType::no_fault:
    return std::make_unique<NoFaultFactory>(drParameters);
  case FrictionLawType::imposed_slip_rates:
    return std::make_unique<ImposedSlipRatesFactory>(drParameters);
  case FrictionLawType::linear_slip_weakening:
    return std::make_unique<LinearSlipWeakeningFactory>(drParameters);
  case FrictionLawType::linear_slip_weakening_forced_rupture_time:
    return std::make_unique<LinearSlipWeakeningForcedRuptureTimeFactory>(drParameters);
  // Prakash-Clifton regularisation for bimaterial faults: see (Pelties et al. 2014)
  case FrictionLawType::linear_slip_weakening_bimaterial:
    return std::make_unique<LinearSlipWeakeningBimaterialFactory>(drParameters);
  case FrictionLawType::rate_and_state_aging_law:
    return std::make_unique<RateAndStateAgingFactory>(drParameters);
  case FrictionLawType::rate_and_state_slip_law:
    return std::make_unique<RateAndStateSlipFactory>(drParameters);
  case FrictionLawType::rate_and_state_velocity_weakening:
    logError() << "friction law 7 currently disabled";
    return std::unique_ptr<AbstractFactory>(nullptr);
  case FrictionLawType::rate_and_state_aging_nucleation:
    logError() << "friction law 101 currently disabled";
    return std::unique_ptr<AbstractFactory>(nullptr);
  case FrictionLawType::rate_and_state_fast_velocity_weakening:
    if (drParameters.isThermalPressureOn == false)
      return std::make_unique<RateAndStateFastVelocityWeakeningFactory>(drParameters);
    else
      return std::make_unique<RateAndStateThermalPressurisationFactory>(drParameters);
  default:
    logError() << "unknown friction law";
    return nullptr;
  }
}

products NoFaultFactory::produce() {
  return {std::make_unique<seissol::initializers::DynamicRupture>(),
          std::make_unique<initializers::NoFaultInitializer>(drParameters),
          std::make_unique<friction_law::NoFault>(drParameters),
          std::make_unique<output::OutputNoFault>(drParameters)};
}

products LinearSlipWeakeningFactory::produce() {
  return {std::make_unique<seissol::initializers::LTS_LinearSlipWeakening>(),
          std::make_unique<initializers::LinearSlipWeakeningInitializer>(drParameters),
          std::make_unique<friction_law::LinearSlipWeakeningLawFL2>(drParameters),
          std::make_unique<output::OutputLinearSlipWeakening>(drParameters)};
}

products RateAndStateAgingFactory::produce() {
  return {std::make_unique<seissol::initializers::LTS_RateAndState>(),
          std::make_unique<initializers::RateAndStateInitializer>(drParameters),
          std::make_unique<friction_law::AgingLaw>(drParameters),
          std::make_unique<output::OutputRateAndState>(drParameters)};
}

products RateAndStateSlipFactory::produce() {
  return {std::make_unique<seissol::initializers::LTS_RateAndState>(),
          std::make_unique<initializers::RateAndStateInitializer>(drParameters),
          std::make_unique<friction_law::SlipLaw>(drParameters),
          std::make_unique<output::OutputRateAndState>(drParameters)};
}

products LinearSlipWeakeningBimaterialFactory::produce() {
  return {std::make_unique<seissol::initializers::LTS_LinearSlipWeakeningBimaterial>(),
          std::make_unique<initializers::LinearSlipWeakeningBimaterialInitializer>(drParameters),
          std::make_unique<friction_law::LinearSlipWeakeningLawBimaterialFL6>(drParameters),
          std::make_unique<output::OutputLinearSlipWeakeningBimaterial>(drParameters)};
}

products LinearSlipWeakeningForcedRuptureTimeFactory::produce() {
  return {
      std::make_unique<seissol::initializers::LTS_LinearSlipWeakeningForcedRuptureTime>(),
      std::make_unique<initializers::LinearSlipWeakeningForcedRuptureTimeInitializer>(drParameters),
      std::make_unique<friction_law::LinearSlipWeakeningLawFL16>(drParameters),
      std::make_unique<output::OutputLinearSlipWeakening>(drParameters)};
}

products ImposedSlipRatesFactory::produce() {
  return {std::make_unique<seissol::initializers::LTS_ImposedSlipRates>(),
          std::make_unique<initializers::ImposedSlipRatesInitializer>(drParameters),
          std::make_unique<friction_law::ImposedSlipRates>(drParameters),
          std::make_unique<output::OutputImposedSlipRates>(drParameters)};
}

products RateAndStateFastVelocityWeakeningFactory::produce() {
  return {std::make_unique<seissol::initializers::LTS_RateAndStateFastVelocityWeakening>(),
          std::make_unique<initializers::RateAndStateFastVelocityInitializer>(drParameters),
          std::make_unique<friction_law::RateAndStateNucFL103>(drParameters),
          std::make_unique<output::OutputRateAndState>(drParameters)};
}

products RateAndStateThermalPressurisationFactory::produce() {
  return {
      std::make_unique<seissol::initializers::LTS_RateAndStateThermalPressurisation>(),
      std::make_unique<initializers::RateAndStateThermalPressurisationInitializer>(drParameters),
      std::make_unique<friction_law::RateAndStateThermalFL103>(drParameters),
      std::make_unique<output::OutputRateAndStateThermalPressurisation>(drParameters)};
}
} // namespace seissol::dr::factory
