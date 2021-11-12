#include "Factory.h"

#include <Solver/Interoperability.h>

namespace seissol::dr::factory {
std::shared_ptr<AbstractFactory> getFactory(dr::DRParameters& drParameters) {
  switch (drParameters.frictionLawType) {
  //! no fault
  case FrictionLawType::no_fault:
    return std::make_shared<NoFaultFactory>(drParameters);
  //! imposed slip rates on the dynamic rupture boundary
  case FrictionLawType::imposed_slip_rates:
    return std::make_shared<ImposedSlipRatesFactory>(drParameters);
  //! coulomb model for linear slip weakening
  case FrictionLawType::linear_slip_weakening:
    return std::make_shared<LinearSlipWeakeningFactory>(drParameters);
  //! linear slip weakening with forced rupture time
  case FrictionLawType::linear_slip_weakening_forced_rupture_time:
    return std::make_shared<LinearSlipWeakeningForcedRuptureTimeFactory>(drParameters);
  //! linear slip weakening with Prakash-Clifton regularisation for bimaterial faults: see (Pelties
  //! et al. 2014)
  case FrictionLawType::linear_slip_weakening_bimaterial:
    return std::make_shared<LinearSlipWeakeningBimaterialFactory>(drParameters);
  //! rate and state aging law
  case FrictionLawType::rate_and_state_aging_law:
    return std::make_shared<RateAndStateAgingFactory>(drParameters);
  //! rate and state slip law
  case FrictionLawType::rate_and_state_slip_law:
    return std::make_shared<RateAndStateSlipFactory>(drParameters);
  //! severe velocity weakening friction as in Ampuero&Ben-Zion2008
  case FrictionLawType::rate_and_state_velocity_weakening:
    logError() << "friction law 7 currently disables";
    return std::shared_ptr<AbstractFactory>(nullptr);
  //! specific conditions for SCEC TPV101 rate and state aging law + nucleation + time and space
  //! dependent
  case FrictionLawType::rate_and_state_aging_nucleation:
    logError() << "friction law 101 currently disables";
    return std::shared_ptr<AbstractFactory>(nullptr);
  //! rate and state fast velocity weakening
  //! with optional thermal pressurisation
  case FrictionLawType::rate_and_state_fast_velocity_weakening:
    if (drParameters.isThermalPressureOn == false)
      //! without thermal pressure
      return std::make_shared<RateAndStateFastVelocityWeakeningFactory>(drParameters);
    else
      //! with thermal pressure (see Noda and Lapusta 2010)
      return std::make_shared<RateAndStateThermalPressurisationFactory>(drParameters);
  default:
    logError() << "unknown friction law";
    return nullptr;
  }
}

products NoFaultFactory::produce() {
  return {std::make_shared<seissol::initializers::DynamicRupture>(),
          std::make_shared<initializers::NoFaultInitializer>(drParameters),
          std::make_shared<friction_law::NoFault>(drParameters),
          std::make_shared<output::OutputNoFault>(drParameters)};
}

products LinearSlipWeakeningFactory::produce() {
  return {std::make_shared<seissol::initializers::LTS_LinearSlipWeakening>(),
          std::make_shared<initializers::LinearSlipWeakeningInitializer>(drParameters),
          std::make_shared<friction_law::LinearSlipWeakeningLawFL2>(drParameters),
          std::make_shared<output::OutputLinearSlipWeakening>(drParameters)};
}

products RateAndStateAgingFactory::produce() {
  return {std::make_shared<seissol::initializers::LTS_RateAndState>(),
          std::make_shared<initializers::RateAndStateInitializer>(drParameters),
          std::make_shared<friction_law::AgingLaw>(drParameters),
          std::make_shared<output::OutputRateAndState>(drParameters)};
}

products RateAndStateSlipFactory::produce() {
  return {std::make_shared<seissol::initializers::LTS_RateAndState>(),
          std::make_shared<initializers::RateAndStateInitializer>(drParameters),
          std::make_shared<friction_law::SlipLaw>(drParameters),
          std::make_shared<output::OutputRateAndState>(drParameters)};
}

products LinearSlipWeakeningBimaterialFactory::produce() {
  return {std::make_shared<seissol::initializers::LTS_LinearSlipWeakeningBimaterial>(),
          std::make_shared<initializers::LinearSlipWeakeningBimaterialInitializer>(drParameters),
          std::make_shared<friction_law::LinearSlipWeakeningLawBimaterialFL6>(drParameters),
          std::make_shared<output::OutputLinearSlipWeakeningBimaterial>(drParameters)};
}

// products Factory_FL_7::produce() {
//  return {new seissol::initializers::LTS_RateAndStateFL3,
//          new initializers::AgingLawInitializer,
//          new friction_law::VelocityWeakening,
//          new output::OutputRateAndStateFL3};
//}

products LinearSlipWeakeningForcedRuptureTimeFactory::produce() {
  return {
      std::make_shared<seissol::initializers::LTS_LinearSlipWeakeningForcedRuptureTime>(),
      std::make_shared<initializers::LinearSlipWeakeningForcedRuptureTimeInitializer>(drParameters),
      std::make_shared<friction_law::LinearSlipWeakeningLawFL16>(drParameters),
      std::make_shared<output::OutputLinearSlipWeakening>(drParameters)};
}

products ImposedSlipRatesFactory::produce() {
  return {std::make_shared<seissol::initializers::LTS_ImposedSlipRates>(),
          std::make_shared<initializers::ImposedSlipRatesInitializer>(drParameters),
          std::make_shared<friction_law::ImposedSlipRates>(drParameters),
          std::make_shared<output::OutputImposedSlipRates>(drParameters)};
}

products RateAndStateFastVelocityWeakeningFactory::produce() {
  return {std::make_shared<seissol::initializers::LTS_RateAndStateFastVelocityWeakening>(),
          std::make_shared<initializers::RateAndStateFastVelocityInitializer>(drParameters),
          std::make_shared<friction_law::RateAndStateNucFL103>(drParameters),
          std::make_shared<output::OutputRateAndState>(drParameters)};
}

products RateAndStateThermalPressurisationFactory::produce() {
  return {
      std::make_shared<seissol::initializers::LTS_RateAndStateThermalPressurisation>(),
      std::make_shared<initializers::RateAndStateThermalPressurisationInitializer>(drParameters),
      std::make_shared<friction_law::RateAndStateThermalFL103>(drParameters),
      std::make_shared<output::OutputRateAndStateThermalPressurisation>(drParameters)};
}
} // namespace seissol::dr::factory
