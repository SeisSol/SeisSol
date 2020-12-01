//
// Created by adrian on 09.07.20.
//

#include <Solver/Interoperability.h>
#include "DR_factory.h"

namespace seissol {
  namespace dr {
    namespace factory {
      AbstractFactory* getFactory(dr::DrParameterT* DynRupParameter) {
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
    }
  }
}
