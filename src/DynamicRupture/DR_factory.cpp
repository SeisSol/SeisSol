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
          case linear_slip_weakening_forced_time_rapture: return new Factory_FL_16;
          //! Coulomb model for linear slip weakening and bimaterial
          case linear_slip_weakening_bimaterial: return new Factory_FL_6;

          //!---- rate and state -----
          //! rate and state aging law
          case rate_and_state_aging_law: return new Factory_FL_3;
          //! rate and state slip law
          case rate_and_state_slip_law: return new Factory_FL_4;
          //! severe velocity weakening friction as in Ampuero&Ben-Zion2008
          case rate_and_state_velocity_weakening: return new Factory_FL_7;
          //! specific conditions for SCEC TPV101 rate and state aging law + time and space dependent nucleation
          case rate_and_state_aging_nucleation:
            throw std::runtime_error("friction law 101 currently disables");
          //! specific conditions for SCEC TPV103 rate and state slip law + time and space dependent nucleation
          case rate_and_state_slip_nucleation:
            if(DynRupParameter->IsTermalPressureOn == false)
              //!without thermal pressure
              return new Factory_FL_103;
            else
              //!with thermal pressure
              return new Factory_FL_103_Thermal;

          default:
              throw std::runtime_error("unknown friction law");
        }
      }
    }
  }
}


int temporary_main() {
    // instantiate all components necessary for the program
    using namespace seissol::dr;
    Friction_law_type FrictionLaw = linear_slip_weakening_forced_time_rapture;

    // in interoperability:
    // seissol::SeisSol::main.getMemoryManager().getDynamicRupture()
    seissol::initializers::LTSTree* dynRupTree = nullptr;
    unsigned * m_ltsFaceToMeshFace;
    std::unordered_map<std::string, double*>  faultParameters;
    seissol::Interoperability interoperability;

    //in memoryManager.h
    //lts::Base* DrLts = nullptr;   //now in seissol::initializers::DynamicRupture implemented
    seissol::initializers::DynamicRupture* DynRup = nullptr;
    seissol::initializers::BaseDrInitializer* DrInitializer = nullptr;
    fr_law::BaseFrictionSolver* FrictonLaw = nullptr;
    output::Output_Base* DrOutput = nullptr;

    //Interoperability:
    //in void seissol::Interoperability::initializeClusteredLts
    //getFL();
    //this allocates m_DrLts, m_DrInitializer, m_FrictonLaw, m_DrOutput
    //seissol::SeisSol::main.getMemoryManager().initializeFrictionFactory(m_FL);


    //in memory manager.cpp
    //void seissol::initializers::MemoryManager::initializeFrictionFactory(Friction_law_type FrictionLaw)
    //factory::AbstractFactory* Factory = factory::getFactory(FrictionLaw);
    //std::tie(DynRup, DrInitializer, FrictonLaw, DrOutput) = Factory->produce();
    //delete Factory;    // prepare the data

    //in memory manager.cpp
    //void seissol::initializers::MemoryManager::fixateLtsTree
    //m_dynRup.addTo(m_dynRupTree);
    DynRup->addTo(*dynRupTree);

    //in interoperability
    //void seissol::Interoperability::initializeCellLocalMatrices()
    DrInitializer->initializeFrictionMatrices(DynRup, dynRupTree, faultParameters, m_ltsFaceToMeshFace, interoperability);

    //computational part:
    //inside TimeCluster.cpp
    //TimeclusterManager.cpp in void seissol::time_stepping::TimeManager::addClusters
    // new TimeCluster with
    //      i_memoryManager.getFrictionLaw(),
    //      i_memoryManager.getDrOutput(),
    // copied to members of TimeCluster:
    //    dr::fr_law::Base* m_FrictonLaw;
    //    dr::output::Base* m_DrOutput,
    // evaluate is called in void seissol::time_stepping::TimeCluster::computeDynamicRupture( seissol::initializers::Layer&  layerData )

    //!FrictonLaw->evaluate(*DynRup /*+ DrLtsTree*/);

    //!DrOutput->tiePointers(*DynRup /*+ DrLtsTree, + faultWriter*/); // pass ptrs of the first cluster    // inside of a compute loop


    //end of program
    DrOutput->postCompute(*DynRup /*+ DrLtsTree*/);    // at the end of program

    //in MemoryManager.h ~MemoryManager()
    delete DynRup;
    delete DrInitializer;
    delete FrictonLaw;
    delete DrOutput;    // DrLtsTree will get deleted as it is done in SeisSol
    return 0;
}