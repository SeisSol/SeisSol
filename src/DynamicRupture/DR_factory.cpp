#include <Solver/Interoperability.h>
#include "DR_factory.h"

seissol::dr::factory::AbstractFactory* seissol::dr::factory::getFactory(dr::DRParameters* DynRupParameter) {
  switch (DynRupParameter->FrictionLawType) {
    //! no fault
    case no_fault: return new Factory_FL_0;
      //! imposed slip rate on the dynamic rupture boundary
    case imposed_slip_rate_on_DR_boundary: return new Factory_FL_33;

      //!---- linear slip weakening -----
      //! coulomb model for linear slip weakening
    case linear_slip_weakening: return new Factory_FL_2;
    case linear_slip_weakening_forced_time_rupture: return new Factory_FL_16;
      //! Coulomb model for linear slip weakening and bimaterial (see pelties 2014)
    case linear_slip_weakening_bimaterial: return new Factory_FL_6;

      //!---- rate and state -----
      //! rate and state aging law
    case rate_and_state_aging_law: return new Factory_FL_3;
      //! rate and state slip law
    case rate_and_state_slip_law: return new Factory_FL_4;
      //! severe velocity weakening friction as in Ampuero&Ben-Zion2008
    case rate_and_state_velocity_weakening: return new Factory_FL_7;
      //! specific conditions for SCEC TPV101 rate and state aging law + nucleation + time and space dependent
    case rate_and_state_aging_nucleation:
      throw std::runtime_error("friction law 101 currently disables");
      //! specific conditions for SCEC TPV103 rate and state slip law + nucleation + time and space dependent parameters
    case rate_and_state_slip_nucleation:
      if(DynRupParameter->IsTermalPressureOn == false)
        //!without thermal pressure
        return new Factory_FL_103;
      else
        //!with thermal pressure (see Noda and Lapusta 2010)
        return new Factory_FL_103_Thermal;

    default:
      throw std::runtime_error("unknown friction law");
  }
}

seissol::dr::factory::products seissol::dr::factory::Factory_FL_0::produce() {
  return {new seissol::initializers::DynamicRupture,
          new seissol::dr::initializers::NoFaultInitializer,
          new seissol::dr::friction_law::NoFault,
          new seissol::dr::output::Output_NoFaultFL0};
};

seissol::dr::factory::products seissol::dr::factory::Factory_FL_2::produce() {
  return {new seissol::initializers::LTS_LinearSlipWeakeningFL2,
          new seissol::dr::initializers::LinearSlipWeakeningFL2Initializer,
          new seissol::dr::friction_law::LinearSlipWeakeningLawFL2,
          new seissol::dr::output::Output_LinearSlipWeakeningFL2};
};

seissol::dr::factory::products seissol::dr::factory::Factory_FL_3::produce() {
  return {new seissol::initializers::LTS_RateAndStateFL3,
          new seissol::dr::initializers::AgingLawInitializer,
          new seissol::dr::friction_law::AgingLaw,
          new seissol::dr::output::Output_RateAndStateFL3};
};

seissol::dr::factory::products seissol::dr::factory::Factory_FL_4::produce() {
  return {new seissol::initializers::LTS_RateAndStateFL3,
          new seissol::dr::initializers::AgingLawInitializer,
          new seissol::dr::friction_law::SlipLaw,
          new seissol::dr::output::Output_RateAndStateFL3};
};

seissol::dr::factory::products seissol::dr::factory::Factory_FL_6::produce() {
  return {new seissol::initializers::LTS_LinearBimaterialFL6,
          new seissol::dr::initializers::LinearBimaterialFL6Initializer,
          new seissol::dr::friction_law::LinearSlipWeakeningLawBimaterialFL6,
          new seissol::dr::output::Output_LinearBimaterialFL6};
};


seissol::dr::factory::products seissol::dr::factory::Factory_FL_7::produce() {
  return {new seissol::initializers::LTS_RateAndStateFL3,
          new seissol::dr::initializers::AgingLawInitializer,
          new seissol::dr::friction_law::VelocityWeakening,
          new seissol::dr::output::Output_RateAndStateFL3};
};

seissol::dr::factory::products seissol::dr::factory::Factory_FL_16::produce() {
  return {new seissol::initializers::LTS_LinearSlipWeakeningFL16,
          new seissol::dr::initializers::LinearSlipWeakeningFL16Initializer,
          new seissol::dr::friction_law::LinearSlipWeakeningLawFL16,
          new seissol::dr::output::Output_LinearSlipWeakeningFL2};
};

seissol::dr::factory::products seissol::dr::factory::Factory_FL_33::produce() {
  return {new seissol::initializers::LTS_ImposedSlipRatesFL33,
          new seissol::dr::initializers::ImposedSlipRatesFL33Initializer,
          new seissol::dr::friction_law::ImposedSlipRates,
          new seissol::dr::output::Output_ImposedSlipRatesFL33};
};

seissol::dr::factory::products seissol::dr::factory::Factory_FL_103::produce() {
  return {new seissol::initializers::LTS_RateAndStateFL103,
          new seissol::dr::initializers::RateAndStateFL103Initializer,
          new seissol::dr::friction_law::RateAndStateNucFL103,
          new seissol::dr::output::Output_RateAndStateFL103};
};

seissol::dr::factory::products seissol::dr::factory::Factory_FL_103_Thermal::produce() {
  return {new seissol::initializers::LTS_RateAndStateFL103TP,
          new seissol::dr::initializers::RateAndStateFL103TPInitializer,
          new seissol::dr::friction_law::RateAndStateThermalFL103,
          new seissol::dr::output::Output_RateAndStateFL103};
};
