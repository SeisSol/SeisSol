//
// Created by adrian on 08.07.20.
//

#ifndef SEISSOL_DR_FACTORY_H
#define SEISSOL_DR_FACTORY_H

#include <iostream>
#include <tuple>
#include <stdexcept>
#include "DR_LTS_Base.h"
#include "DR_initializer_base.h"
#include "DR_friction_law.h"
#include "DR_output.h"


namespace seissol {
    namespace dr {
        namespace factory {    using products = std::tuple<lts::Base*, initializer::Base*, fr_law::Base*, output::Base*>;    class AbstractFactory {
            public:
                virtual ~AbstractFactory() {}
                virtual products produce() = 0;
            };    class FL_16 : public AbstractFactory {
                virtual products produce() override {
                    return std::make_tuple(new lts::FL_16,
                                           new initializer::FL_16,
                                           new fr_law::FL_16,
                                           new output::FL_16);
                }
            };    class FL_17 : public AbstractFactory {
                virtual products produce() override {
                    return std::make_tuple(new lts::FL_16,
                                           new initializer::FL_16,
                                           new fr_law::FL_17,
                                           new output::FL_16);
                }
            };    class FL_33 : public AbstractFactory {
                virtual products produce() override {
                    return std::make_tuple(new lts::FL_33,
                                           new initializer::FL_33,
                                           new fr_law::FL_33,
                                           new output::FL_33);
                }
            };    AbstractFactory* getFactory(Friction_law_type FrictionLawID) {
                switch (FrictionLawID) {
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

int main() {
    // instantiate all components necessary for the program
    using namespace seissol::dr;
    Friction_law_type FrictionLaw = Linear_slip_weakening_forced_time_rapture;

    // in interoperability:
    // seissol::SeisSol::main.getMemoryManager().getDynamicRupture()
    seissol::initializers::LTSTree* dynRupTree = nullptr;

    //in memoryManager.h
    lts::Base* DrLts = nullptr;
    initializer::Base* DrInitializer = nullptr;
    fr_law::Base* FrictonLaw = nullptr;
    output::Base* DrOutput = nullptr;

    //Interoperability:
    //in void seissol::Interoperability::initializeClusteredLts
    //getFL();
    //seissol::SeisSol::main.getMemoryManager().initializeFrictionFactory(m_FL);


    //in memory manager.cpp
    //void seissol::initializers::MemoryManager::initializeFrictionFactory(Friction_law_type FrictionLaw)
    factory::AbstractFactory* Factory = factory::getFactory(FrictionLaw);
    std::tie(DrLts, DrInitializer, FrictonLaw, DrOutput) = Factory->produce();
    delete Factory;    // prepare the data

    //in memory manager.cpp
    //void seissol::initializers::MemoryManager::fixateLtsTree
    DrLts->addVars(*dynRupTree /*+ DrLtsTree*/);
    DrInitializer->initializeFrictionMatrices(DrLts, dynRupTree /*+ DrLtsTree, + something from Easy*/);

    //TODO: Dont know where to put this
    //computational part:
    DrOutput->tiePointers(*DrLts /*+ DrLtsTree, + faultWriter*/); // pass ptrs of the first cluster    // inside of a compute loop


    //inside TimeCluster.cpp
    //TimeclusterManager.cpp in void seissol::time_stepping::TimeManager::addClusters
    // new TimeCluster with
    //      i_memoryManager.getFrictionLaw(),
    //      i_memoryManager.getDrLts(),
    // copied to members of TimeCluster:
    //    dr::fr_law::Base* m_FrictonLaw;
    //    dr::lts::Base* m_DrLts;
    // evaluate is called in void seissol::time_stepping::TimeCluster::computeDynamicRupture( seissol::initializers::Layer&  layerData )
    FrictonLaw->evaluate(*DrLts /*+ DrLtsTree*/);

    //TODO: dont know where to put these
    //end of program
    DrOutput->postCompute(*DrLts /*+ DrLtsTree*/);    // at the end of program
    delete DrLts;
    delete DrInitializer;
    delete FrictonLaw;
    delete DrOutput;    // DrLtsTree will get deleted as it is done in SeisSol
    return 0;
}













#endif //SEISSOL_DR_FACTORY_H
