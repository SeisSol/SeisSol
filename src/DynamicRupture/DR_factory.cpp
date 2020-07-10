//
// Created by adrian on 09.07.20.
//

#include "DR_factory.h"

namespace seissol {
    namespace dr {
        namespace factory {
            AbstractFactory* getFactory(Friction_law_type FrictionLawID) {
                switch (FrictionLawID) {
                    case Linear_slip_weakening: return new FL_2;
                    case Linear_slip_weakening_forced_time_rapture: return new FL_16;
                    case 17: return new FL_17;
                    case 33: return new FL_33;
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
    Friction_law_type FrictionLaw = Linear_slip_weakening_forced_time_rapture;

    // in interoperability:
    // seissol::SeisSol::main.getMemoryManager().getDynamicRupture()
    seissol::initializers::LTSTree* dynRupTree = nullptr;

    //in memoryManager.h
    //lts::Base* DrLts = nullptr;   //now in seissol::initializers::DynamicRupture implemented
    seissol::initializers::DynamicRupture* DynRup = nullptr;
    initializer::Base* DrInitializer = nullptr;
    fr_law::Base* FrictonLaw = nullptr;
    output::Base* DrOutput = nullptr;

    //Interoperability:
    //in void seissol::Interoperability::initializeClusteredLts
    //getFL();
    //this allocates m_DrLts, m_DrInitializer, m_FrictonLaw, m_DrOutput
    //seissol::SeisSol::main.getMemoryManager().initializeFrictionFactory(m_FL);


    //in memory manager.cpp
    //void seissol::initializers::MemoryManager::initializeFrictionFactory(Friction_law_type FrictionLaw)
    factory::AbstractFactory* Factory = factory::getFactory(FrictionLaw);
    std::tie(DynRup, DrInitializer, FrictonLaw, DrOutput) = Factory->produce();
    delete Factory;    // prepare the data

    //in memory manager.cpp
    //void seissol::initializers::MemoryManager::fixateLtsTree
    //m_dynRup.addTo(m_dynRupTree);
    DynRup->addTo(*dynRupTree);

    //in interoperability
    //void seissol::Interoperability::initializeCellLocalMatrices()
    DrInitializer->initializeFrictionMatrices(DynRup, dynRupTree /*+ DrLtsTree, + something from Easy*/);

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
    FrictonLaw->evaluate(*DynRup /*+ DrLtsTree*/);
    DrOutput->tiePointers(*DynRup /*+ DrLtsTree, + faultWriter*/); // pass ptrs of the first cluster    // inside of a compute loop


    //end of program
    //TODO: dont know where to put this
    DrOutput->postCompute(*DynRup /*+ DrLtsTree*/);    // at the end of program

    //in MemoryManager.h ~MemoryManager()
    delete DynRup;
    delete DrInitializer;
    delete FrictonLaw;
    delete DrOutput;    // DrLtsTree will get deleted as it is done in SeisSol
    return 0;
}