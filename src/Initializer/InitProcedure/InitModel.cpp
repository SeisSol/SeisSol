
#include <vector>
#include "Initializer/ParameterDB.h"
#include "Initializer/InputParameters.hpp"
#include "Initializer/CellLocalMatrices.h"
#include "Initializer/LTS.h"
#include "Initializer/tree/LTSTree.hpp"
#include "Initializer/time_stepping/common.hpp"
#include "Physics/Attenuation.hpp"
#include "Equations/datastructures.hpp"
#include "Initializer/tree/Lut.hpp"
#include "Initializer/typedefs.hpp"

#include "SeisSol.h"
#include "Init.hpp"
#include "InitModel.hpp"

#include "Parallel/MPI.h"

#include <cmath>

using namespace seissol::initializer;

using MaterialClass = seissol::model::MaterialClass;
using Plasticity = seissol::model::Plasticity;

template<typename T>
void synchronize(const seissol::initializers::Variable<T>& handle)
{
  initializers::MemoryManager& memoryManager = seissol::SeisSol::main.getMemoryManager();
  unsigned *const (&meshToLts)[seissol::initializers::Lut::MaxDuplicates] = memoryManager.getLtsLut()->getMeshToLtsLut(handle.mask);
  unsigned* duplicatedMeshIds = memoryManager.getLtsLut()->getDuplicatedMeshIds(handle.mask);
  unsigned numberOfDuplicatedMeshIds = memoryManager.getLtsLut()->getNumberOfDuplicatedMeshIds(handle.mask);
  T* var = memoryManager.getLtsTree()->var(handle);
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned dupMeshId = 0; dupMeshId < numberOfDuplicatedMeshIds; ++dupMeshId) {
    unsigned meshId = duplicatedMeshIds[dupMeshId];
    T* ref = &var[ meshToLts[0][meshId] ];
    for (unsigned dup = 1; dup < seissol::initializers::Lut::MaxDuplicates && meshToLts[dup][meshId] != std::numeric_limits<unsigned>::max(); ++dup) {
      memcpy(reinterpret_cast<void*>(&var[ meshToLts[dup][meshId] ]),
             reinterpret_cast<void*>(ref),
             sizeof(T));
    }
  }
}

template<typename T>
std::vector<T> queryDB(seissol::initializers::QueryGenerator* queryGen, const std::string& fileName, size_t size) {
    std::vector<T> vectorDB(size);
    seissol::initializers::MaterialParameterDB<T> parameterDB;
    parameterDB.setMaterialVector(&vectorDB);
    parameterDB.evaluateModel(fileName, queryGen);
    return vectorDB;
}

template<typename T>
void copyInit(T& target, const T& value) {
  new (&target) T (value);
}

void initializeCellMaterial() {
    const auto& ssp = seissol::SeisSol::main.getSeisSolParameters();
    const auto& meshReader = seissol::SeisSol::main.meshReader();
    initializers::MemoryManager& memoryManager = seissol::SeisSol::main.getMemoryManager();

    // unpack ghost layer (merely a re-ordering operation)
    std::vector<std::array<std::array<double, 3>, 4>> ghostVertices;
    std::vector<int> ghostMaterials;
    std::unordered_map<int, std::vector<unsigned>> ghostIdxMap;
    for (const auto& neighbor : meshReader.getMPINeighborVertices()) {
      ghostIdxMap[neighbor.first].reserve(neighbor.second.size());
      for (const auto& vertices : neighbor.second) {
        ghostIdxMap[neighbor.first].push_back(ghostVertices.size());
        ghostVertices.push_back(vertices);
        ghostMaterials.push_back(0); // TODO: change that
      }
    }

    seissol::initializers::QueryGenerator* queryGen = seissol::initializers::getBestQueryGenerator(
      seissol::initializer::parameters::modelAnelastic(),
      ssp.model.plasticity,
      seissol::initializer::parameters::modelAnisotropic(),
      seissol::initializer::parameters::modelPoroelastic(),
      ssp.model.useCellHomogenizedMaterial,
      seissol::initializers::C2VArray::fromMeshReader(meshReader));
    seissol::initializers::QueryGenerator* queryGenGhost = seissol::initializers::getBestQueryGenerator(
      seissol::initializer::parameters::modelAnelastic(),
      ssp.model.plasticity,
      seissol::initializer::parameters::modelAnisotropic(),
      seissol::initializer::parameters::modelPoroelastic(),
      ssp.model.useCellHomogenizedMaterial,
      seissol::initializers::C2VArray::fromVectors(ghostVertices, ghostMaterials));
    auto materialsDB = queryDB<MaterialClass>(queryGen, ssp.model.materialFileName, meshReader.getElements().size());
    auto materialsDBGhost = queryDB<MaterialClass>(queryGenGhost, ssp.model.materialFileName, ghostVertices.size());
    std::vector<Plasticity> plasticityDB;
    if (ssp.model.plasticity) {
      // plasticity information is only needed on all interior+copy cells.
        plasticityDB = queryDB<Plasticity>(queryGen, ssp.model.materialFileName, meshReader.getElements().size());
    }

#if defined(USE_VISCOELASTIC) || defined(USE_VISCOELASTIC2)
    // we need to compute all model parameters before we can use them...
    // TODO: integrate this with the Viscoelastic material class or the ParameterDB directly?
    logDebug() << "Initializing attenuation.";
#ifdef OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < materialsDB.size(); ++i) {
        auto& cellmat = materialsDB[i];
        seissol::physics::fitAttenuation(cellmat, ssp.model.freqCentral, ssp.model.freqRatio);
    }
#ifdef OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < materialsDBGhost.size(); ++i) {
        auto& cellmat = materialsDBGhost[i];
        seissol::physics::fitAttenuation(cellmat, ssp.model.freqCentral, ssp.model.freqRatio);
    }
#endif

  logDebug() << "Setting cell materials in the LTS tree (for interior and copy layers).";
  const auto& elements = meshReader.getElements();
  unsigned* ltsToMesh = memoryManager.getLtsLut()->getLtsToMeshLut(memoryManager.getLts()->material.mask);

  for (seissol::initializers::LTSTree::leaf_iterator it = memoryManager.getLtsTree()->beginLeaf(seissol::initializers::LayerMask(Ghost)); it != memoryManager.getLtsTree()->endLeaf(); ++it) {
    auto* cellInformation = it->var(memoryManager.getLts()->cellInformation);
    auto* materialArray = it->var(memoryManager.getLts()->material);
    auto* plasticityArray = ssp.model.plasticity ? it->var(memoryManager.getLts()->plasticity) : nullptr;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t cell = 0; cell < it->getNumberOfCells(); ++cell) {
      // set the materials for the cell volume and its faces
      auto meshId = ltsToMesh[cell];
      auto& material = materialArray[cell];
      const auto& localMaterial = materialsDB[meshId];
      const auto& element = elements[meshId];
      const auto& localCellInformation = cellInformation[cell];

      copyInit(material.local, localMaterial);
      for (std::size_t side = 0; side < 4; ++side) {
          if (isInteriorFaceType(localCellInformation.faceTypes[side])) {
              // use the neighbor face material info in case that we are not at a boundary
              if (element.neighborRanks[side] == seissol::MPI::mpi.rank()) {
                // material from interior or copy
                auto neighbor = element.neighbors[side];
                copyInit(material.neighbor[side], materialsDB[neighbor]);
              }
              else {
                // material from ghost layer (computed locally)
                auto neighborRank = element.neighborRanks[side];
                auto neighborRankIdx = element.mpiIndices[side];
                auto materialGhostIdx = ghostIdxMap.at(neighborRank)[neighborRankIdx];
                copyInit(material.neighbor[side], materialsDBGhost[materialGhostIdx]);
              }
          }
          else {
              // otherwise, use the material from the own cell (TODO: is this really the best thing for boundary conditions?)
              copyInit(material.neighbor[side], localMaterial);
          }
        }
      
      // if enabled, set up the plasticity as well
      // TODO: move to material initalization maybe? Or an initializer for the PlasticityData struct?
        if (ssp.model.plasticity) {
          auto& plasticity = plasticityArray[cell];
            const auto& localPlasticity = plasticityDB[meshId];

            plasticity.initialLoading[0] = localPlasticity.s_xx;
            plasticity.initialLoading[1] = localPlasticity.s_yy;
            plasticity.initialLoading[2] = localPlasticity.s_zz;
            plasticity.initialLoading[3] = localPlasticity.s_xy;
            plasticity.initialLoading[4] = localPlasticity.s_yz;
            plasticity.initialLoading[5] = localPlasticity.s_xz;

            double angularFriction = std::atan(localPlasticity.bulkFriction);

            plasticity.cohesionTimesCosAngularFriction = localPlasticity.plastCo * std::cos(angularFriction);
            plasticity.sinAngularFriction = std::sin(angularFriction);
#ifndef USE_ANISOTROPIC
            plasticity.mufactor = 1.0 / (2.0 * material.local.mu);
#else
            plasticity.mufactor = 3.0 / (2.0 * (material.local.c44 + material.local.c55 + material.local.c66));
#endif
        }
    }
    ltsToMesh += it->getNumberOfCells();
  }

    // set tv for all time clusters
    if (ssp.model.plasticity) {
      seissol::SeisSol::main.timeManager().setTv(ssp.model.tv);
    }

    // TODO: check if sync is still necessary
    // synchronize data
    synchronize(memoryManager.getLts()->material);
    if (ssp.model.plasticity) {
      synchronize(memoryManager.getLts()->plasticity);
    }
}

struct LtsInfo {
    unsigned* ltsMeshToFace = nullptr;
    MeshStructure* meshStructure = nullptr;
    TimeStepping timeStepping;

    // IMPORTANT: DO NOT DEALLOCATE THE ABOVE POINTERS... THEY ARE PASSED ON AND REQUIRED DURING RUNTIME
};

void initializeCellMatrices(LtsInfo& ltsInfo)
{
  const auto& ssp = seissol::SeisSol::main.getSeisSolParameters();

  // \todo Move this to some common initialization place
  auto& meshReader = seissol::SeisSol::main.meshReader();
  auto& memoryManager = seissol::SeisSol::main.getMemoryManager();

  seissol::initializers::initializeCellLocalMatrices( meshReader,
                                                      memoryManager.getLtsTree(),
                                                      memoryManager.getLts(),
                                                      memoryManager.getLtsLut(),
                                                      ltsInfo.timeStepping);

  seissol::initializers::initializeDynamicRuptureMatrices( meshReader,
                                                           memoryManager.getLtsTree(),
                                                           memoryManager.getLts(),
                                                           memoryManager.getLtsLut(),
                                                           memoryManager.getDynamicRuptureTree(),
                                                           memoryManager.getDynamicRupture(),
                                                           ltsInfo.ltsMeshToFace,
                                                           *memoryManager.getGlobalDataOnHost(),
                                                           ltsInfo.timeStepping );

  memoryManager.initFrictionData();
  // seissol::SeisSol::main.getMemoryManager().getFaultOutputManager()->initFaceToLtsMap(); // moved

  seissol::initializers::initializeBoundaryMappings(meshReader,
                                                    memoryManager.getEasiBoundaryReader(),
                                                    memoryManager.getLtsTree(),
                                                    memoryManager.getLts(),
                                                    memoryManager.getLtsLut());

#ifdef ACL_DEVICE
  initializers::copyCellMatricesToDevice(memoryManager.getLtsTree(),
                                         memoryManager.getLts(),
                                         memoryManager.getDynamicRuptureTree(),
                                         memoryManager.getDynamicRupture(),
                                         memoryManager.getBoundaryTree(),
                                         memoryManager.getBoundary());

  memoryManager.recordExecutionPaths(ssp.model.plasticity);
#endif

  // synchronize data
  synchronize(memoryManager.getLts()->dofs);
  if (kernels::size<tensor::Qane>() > 0) {
    synchronize(memoryManager.getLts()->dofsAne);
  }
}

void initializeClusteredLts(LtsInfo& ltsInfo) {
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

  // derive lts setups
  seissol::initializers::time_stepping::deriveLtsSetups( ltsInfo.timeStepping.numberOfLocalClusters,
                                                         ltsInfo.meshStructure,
                                                         m_ltsTree->var(m_lts->cellInformation) );

}

void initializeMemoryLayout(LtsInfo& ltsInfo) {
  const auto& ssp = seissol::SeisSol::main.getSeisSolParameters();

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


void seissol::initializer::initprocedure::initModel() {
    SCOREP_USER_REGION("init_model", SCOREP_USER_REGION_TYPE_FUNCTION);

    logInfo(seissol::MPI::mpi.rank()) << "Begin init model.";

    // Call the pre mesh initialization hook
	seissol::Modules::callHook<seissol::PRE_MODEL>();

    seissol::Stopwatch watch;
	watch.start();

  LtsInfo ltsInfo;
  
  // these four methods need to be called in this order.

  // init LTS
  logInfo(seissol::MPI::mpi.rank()) << "Initialize LTS.";
  initializeClusteredLts(ltsInfo);

  // init cell materials (needs LTS, to place the material in; this part was translated from FORTRAN)
  logInfo(seissol::MPI::mpi.rank()) << "Initialize cell material parameters.";
	initializeCellMaterial();

  // init memory layout (needs cell material values to initialize e.g. displacements correctly)
  logInfo(seissol::MPI::mpi.rank()) << "Initialize Memory layout.";
  initializeMemoryLayout(ltsInfo);

  // init cell matrices (TODO: what is the dependency here? However, it was called after init memory layout in the FORTRAN code, so we'll do this here as well)
  logInfo(seissol::MPI::mpi.rank()) << "Initialize cell-local matrices.";
  initializeCellMatrices(ltsInfo);

	watch.pause();
	watch.printTime("Model initialized in:");

    // Call the post mesh initialization hook
	seissol::Modules::callHook<seissol::POST_MODEL>();

  logInfo(seissol::MPI::mpi.rank()) << "End init model.";
}
