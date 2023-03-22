#include "SeisSol.h"
#include "Init.hpp"
#include "InitLts.hpp"

#include "Initializer/time_stepping/common.hpp"
#include "Initializer/time_stepping/LtsLayout.h"
#include "Initializer/MemoryManager.h"

void initializeClusteredLts(seissol::initializer::initprocedure::LtsInfo& ltsInfo) {
  const auto& ssp = seissol::SeisSol::main.getSeisSolParameters();

  // assert a valid clustering
  assert(ssp.timestepping.lts.rate > 0 );

  // either derive a GTS or LTS layout
  if(ssp.timestepping.lts.rate == 1 ) {
    seissol::SeisSol::main.getLtsLayout().deriveLayout( single, 1);
  }
  else {
    seissol::SeisSol::main.getLtsLayout().deriveLayout(multiRate, ssp.timestepping.lts.rate );
  }

  // get the mesh structure
  seissol::SeisSol::main.getLtsLayout().getMeshStructure( ltsInfo.meshStructure );

  // get time stepping
  seissol::SeisSol::main.getLtsLayout().getCrossClusterTimeStepping( ltsInfo.timeStepping );

  unsigned* numberOfDRCopyFaces;
  unsigned* numberOfDRInteriorFaces;
  // get cell information & mappings
  seissol::SeisSol::main.getLtsLayout().getDynamicRuptureInformation( ltsInfo.ltsMeshToFace,
                                                                      numberOfDRCopyFaces,
                                                                      numberOfDRInteriorFaces );

  // swapped with the previous method, but that should have no effect
  seissol::SeisSol::main.getMemoryManager().initializeFrictionLaw();

  seissol::SeisSol::main.getMemoryManager().fixateLtsTree(ltsInfo.timeStepping,
                                                          ltsInfo.meshStructure,
                                                          numberOfDRCopyFaces,
                                                          numberOfDRInteriorFaces,
                                                          ssp.model.plasticity);

  delete[] numberOfDRCopyFaces;
  delete[] numberOfDRInteriorFaces;

  const auto& m_ltsTree = seissol::SeisSol::main.getMemoryManager().getLtsTree();
  const auto& m_lts = seissol::SeisSol::main.getMemoryManager().getLts();

  unsigned* ltsToMesh;
  unsigned numberOfMeshCells;
  // get cell information & mappings
  seissol::SeisSol::main.getLtsLayout().getCellInformation( m_ltsTree->var(m_lts->cellInformation),
                                                            ltsToMesh,
                                                            numberOfMeshCells );

  // TODO: move all of this method to the MemoryManager
  seissol::SeisSol::main.getMemoryManager().getLtsLutUnsafe().createLuts(  m_ltsTree,
                        ltsToMesh,
                        numberOfMeshCells );

  delete[] ltsToMesh;

  const auto& m_dynRupTree = seissol::SeisSol::main.getMemoryManager().getDynamicRuptureTree();
  const auto& m_dynRup = seissol::SeisSol::main.getMemoryManager().getDynamicRupture();

  // derive lts setups
  seissol::initializers::time_stepping::deriveLtsSetups( ltsInfo.timeStepping.numberOfLocalClusters,
                                                         ltsInfo.meshStructure,
                                                         m_ltsTree->var(m_lts->cellInformation) );

  // init memory layout (method break before here; including setting materials first... It should not cause any problems if we merge these methods... I hope.)

  // initialize memory layout
  seissol::SeisSol::main.getMemoryManager().initializeMemoryLayout();

  // add clusters
  seissol::SeisSol::main.timeManager().addClusters(ltsInfo.timeStepping,
                                                   ltsInfo.meshStructure,
                                                   seissol::SeisSol::main.getMemoryManager(),
                                                   ssp.model.plasticity);

  // get backward coupling
  // m_globalData = seissol::SeisSol::main.getMemoryManager().getGlobalDataOnHost();


  // initialize face lts trees
  seissol::SeisSol::main.getMemoryManager().fixateBoundaryLtsTree();
}

void seissol::initializer::initprocedure::initLts(seissol::initializer::initprocedure::LtsInfo& ltsInfo) {
    SCOREP_USER_REGION("init_lts", SCOREP_USER_REGION_TYPE_FUNCTION);

    logInfo() << "Begin init LTS.";

    // Call the pre mesh initialization hook
	seissol::Modules::callHook<seissol::PRE_MODEL>();

    seissol::Stopwatch watch;
	watch.start();

    logInfo() << "Initialize LTS.";
	initializeClusteredLts(ltsInfo);

	watch.pause();
	watch.printTime("LTS initialized in:");

    // Call the post mesh initialization hook
	seissol::Modules::callHook<seissol::POST_MODEL>();

    logInfo() << "End init LTS.";
}
