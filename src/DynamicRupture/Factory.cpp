#include "Factory.h"

#include <Solver/Interoperability.h>

namespace seissol::dr::factory {
AbstractFactory* getFactory(dr::DRParameters* DynRupParameter) {
  switch (DynRupParameter->frictionLawType) {
  //! no fault
  case FrictionLawType::no_fault:
    return new NoFaultFactory;
  //! imposed slip rates on the dynamic rupture boundary
  case FrictionLawType::imposed_slip_rates:
    return new ImposedSlipRatesFactory;
  //! coulomb model for linear slip weakening
  case FrictionLawType::linear_slip_weakening:
    return new LinearSlipWeakeningFactory;
  //! linear slip weakening with forced rupture time
  case FrictionLawType::linear_slip_weakening_forced_rupture_time:
    return new LinearSlipWeakeningForcedRuptureTimeFactory;
  //! linear slip weakening with Prakash-Clifton regularisation for bimaterial faults: see (Pelties
  //! et al. 2014)
  case FrictionLawType::linear_slip_weakening_bimaterial:
    return new LinearSlipWeakeningBimaterialFactory;
  //! rate and state aging law
  case FrictionLawType::rate_and_state_aging_law:
    return new RateAndStateAgingFactory;
  //! rate and state slip law
  case FrictionLawType::rate_and_state_slip_law:
    return new RateAndStateSlipFactory;
  //! severe velocity weakening friction as in Ampuero&Ben-Zion2008
  case FrictionLawType::rate_and_state_velocity_weakening:
    throw std::runtime_error("friction law 7 currently disables");
    // return new Factory_FL_7;
  //! specific conditions for SCEC TPV101 rate and state aging law + nucleation + time and space
  //! dependent
  case FrictionLawType::rate_and_state_aging_nucleation:
    throw std::runtime_error("friction law 101 currently disables");
  //! rate and state fast velocity weakening
  //! with optional thermal pressurisation
  case FrictionLawType::rate_and_state_fast_velocity_weakening:
    if (DynRupParameter->isThermalPressureOn == false)
      //! without thermal pressure
      return new RateAndStateFastVelocityWeakeningFactory;
    else
      //! with thermal pressure (see Noda and Lapusta 2010)
      return new RateAndStateThermalPressurisationFactory;
  default:
    throw std::runtime_error("unknown friction law");
  }
}

products NoFaultFactory::produce() {
  return {new seissol::initializers::DynamicRupture,
          new initializers::NoFaultInitializer,
          new friction_law::NoFault,
          new output::OutputNoFault};
}

products LinearSlipWeakeningFactory::produce() {
  return {new seissol::initializers::LTS_LinearSlipWeakening,
          new initializers::LinearSlipWeakeningFL2Initializer,
          new friction_law::LinearSlipWeakeningLawFL2,
          new output::OutputLinearSlipWeakening};
}

products RateAndStateAgingFactory::produce() {
  return {new seissol::initializers::LTS_RateAndState,
          new initializers::AgingLawInitializer,
          new friction_law::AgingLaw,
          new output::OutputRateAndState};
}

products RateAndStateSlipFactory::produce() {
  return {new seissol::initializers::LTS_RateAndState,
          new initializers::AgingLawInitializer,
          new friction_law::SlipLaw,
          new output::OutputRateAndState};
}

products LinearSlipWeakeningBimaterialFactory::produce() {
  return {new seissol::initializers::LTS_LinearSlipWeakeningBimaterial,
          new initializers::LinearBimaterialFL6Initializer,
          new friction_law::LinearSlipWeakeningLawBimaterialFL6,
          new output::OutputLinearSlipWeakeningBimaterial};
}

// products Factory_FL_7::produce() {
//  return {new seissol::initializers::LTS_RateAndStateFL3,
//          new initializers::AgingLawInitializer,
//          new friction_law::VelocityWeakening,
//          new output::OutputRateAndStateFL3};
//}

products LinearSlipWeakeningForcedRuptureTimeFactory::produce() {
  return {new seissol::initializers::LTS_LinearSlipWeakeningForcedRuptureTime,
          new initializers::LinearSlipWeakeningFL16Initializer,
          new friction_law::LinearSlipWeakeningLawFL16,
          new output::OutputLinearSlipWeakening};
}

products ImposedSlipRatesFactory::produce() {
  return {new seissol::initializers::LTS_ImposedSlipRates,
          new initializers::ImposedSlipRatesFL33Initializer,
          new friction_law::ImposedSlipRates,
          new output::OutputImposedSlipRates};
}

products RateAndStateFastVelocityWeakeningFactory::produce() {
  return {new seissol::initializers::LTS_RateAndStateFastVelocityWeakening,
          new initializers::RateAndStateFL103Initializer,
          new friction_law::RateAndStateNucFL103,
          new output::OutputRateAndStateFastVelocityWeakening};
}

products RateAndStateThermalPressurisationFactory::produce() {
  return {new seissol::initializers::LTS_RateAndStateThermalPressurisation,
          new initializers::RateAndStateFL103TPInitializer,
          new friction_law::RateAndStateThermalFL103,
          new output::OutputRateAndStateFastVelocityWeakening};
}
} // namespace seissol::dr::factory
