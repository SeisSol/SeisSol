
#include "Equations/Datastructures.h"
#include "Initializer/CellLocalMatrices.h"
#include "Initializer/LTS.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/TimeStepping/Common.h"
#include "Initializer/Tree/LTSSync.h"
#include "Initializer/Tree/LTSTree.h"
#include "Initializer/Tree/Lut.h"
#include "Initializer/Typedefs.h"
#include "Physics/Attenuation.h"
#include <Initializer/BasicTypedefs.h>
#include <Initializer/MemoryManager.h>
#include <Initializer/Parameters/ModelParameters.h>
#include <Initializer/Tree/Layer.h>
#include <Model/CommonDatastructures.h>
#include <Model/Plasticity.h>
#include <Modules/Modules.h>
#include <Monitoring/Stopwatch.h>
#include <Physics/InstantaneousTimeMirrorManager.h>
#include <array>
#include <cassert>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <utils/logger.h>
#include <vector>

#include "InitModel.h"
#include "SeisSol.h"

#include "Parallel/MPI.h"

#include <cmath>

using namespace seissol::initializer;

namespace {

using MaterialT = seissol::model::MaterialT;
using Plasticity = seissol::model::Plasticity;

template <typename T>
std::vector<T> queryDB(seissol::initializer::QueryGenerator* queryGen,
                       const std::string& fileName,
                       size_t size) {
  std::vector<T> vectorDB(size);
  seissol::initializer::MaterialParameterDB<T> parameterDB;
  parameterDB.setMaterialVector(&vectorDB);
  parameterDB.evaluateModel(fileName, queryGen);
  return vectorDB;
}

void initializeCellMaterial(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  const auto& meshReader = seissolInstance.meshReader();
  initializer::MemoryManager& memoryManager = seissolInstance.getMemoryManager();

  // unpack ghost layer (merely a re-ordering operation, since the CellToVertexArray right now
  // requires an vector there)
  std::vector<std::array<std::array<double, 3>, 4>> ghostVertices;
  std::vector<int> ghostGroups;
  std::unordered_map<int, std::vector<unsigned>> ghostIdxMap;
  for (const auto& neighbor : meshReader.getGhostlayerMetadata()) {
    ghostIdxMap[neighbor.first].reserve(neighbor.second.size());
    for (const auto& metadata : neighbor.second) {
      ghostIdxMap[neighbor.first].push_back(ghostVertices.size());
      std::array<std::array<double, 3>, 4> vertices{};
      for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          vertices[i][j] = metadata.vertices[i][j];
        }
      }
      ghostVertices.emplace_back(vertices);
      ghostGroups.push_back(metadata.group);
    }
  }

  // just a helper function for better readability
  auto getBestQueryGenerator = [&](const seissol::initializer::CellToVertexArray& ctvArray) {
    return seissol::initializer::getBestQueryGenerator(
        seissol::initializer::parameters::isModelAnelastic(),
        seissolParams.model.plasticity,
        seissol::initializer::parameters::isModelAnisotropic(),
        seissol::initializer::parameters::isModelPoroelastic(),
        seissolParams.model.useCellHomogenizedMaterial,
        ctvArray);
  };

  // material retrieval for copy+interior layers
  seissol::initializer::QueryGenerator* queryGen =
      getBestQueryGenerator(seissol::initializer::CellToVertexArray::fromMeshReader(meshReader));
  auto materialsDB = queryDB<MaterialT>(
      queryGen, seissolParams.model.materialFileName, meshReader.getElements().size());

  // plasticity (if needed)
  std::vector<Plasticity> plasticityDB;
  if (seissolParams.model.plasticity) {
    // plasticity information is only needed on all interior+copy cells.
    plasticityDB = queryDB<Plasticity>(
        queryGen, seissolParams.model.materialFileName, meshReader.getElements().size());
  }

  // material retrieval for ghost layers
  seissol::initializer::QueryGenerator* queryGenGhost = getBestQueryGenerator(
      seissol::initializer::CellToVertexArray::fromVectors(ghostVertices, ghostGroups));
  auto materialsDBGhost =
      queryDB<MaterialT>(queryGenGhost, seissolParams.model.materialFileName, ghostVertices.size());

#if defined(USE_VISCOELASTIC) || defined(USE_VISCOELASTIC2)
  // we need to compute all model parameters before we can use them...
  // TODO(David): integrate this with the Viscoelastic material class or the ParameterDB directly?
  logDebug() << "Initializing attenuation.";
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (size_t i = 0; i < materialsDB.size(); ++i) {
    auto& cellMat = materialsDB[i];
    seissol::physics::fitAttenuation(
        cellMat, seissolParams.model.freqCentral, seissolParams.model.freqRatio);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (size_t i = 0; i < materialsDBGhost.size(); ++i) {
    auto& cellMat = materialsDBGhost[i];
    seissol::physics::fitAttenuation(
        cellMat, seissolParams.model.freqCentral, seissolParams.model.freqRatio);
  }
#endif

  logDebug() << "Setting cell materials in the LTS tree (for interior and copy layers).";
  const auto& elements = meshReader.getElements();
  const unsigned* ltsToMesh =
      memoryManager.getLtsLut()->getLtsToMeshLut(memoryManager.getLts()->material.mask);

  for (auto& layer : memoryManager.getLtsTree()->leaves(Ghost)) {
    auto* cellInformation = layer.var(memoryManager.getLts()->cellInformation);
    auto* materialArray = layer.var(memoryManager.getLts()->material);
    auto* plasticityArray =
        seissolParams.model.plasticity ? layer.var(memoryManager.getLts()->plasticity) : nullptr;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t cell = 0; cell < layer.getNumberOfCells(); ++cell) {
      // set the materials for the cell volume and its faces
      auto meshId = ltsToMesh[cell];
      auto& material = materialArray[cell];
      const auto& localMaterial = materialsDB[meshId];
      const auto& element = elements[meshId];
      const auto& localCellInformation = cellInformation[cell];

      initAssign(material.local, localMaterial);
      for (std::size_t side = 0; side < 4; ++side) {
        if (isInternalFaceType(localCellInformation.faceTypes[side])) {
          // use the neighbor face material info in case that we are not at a boundary
          if (element.neighborRanks[side] == seissol::MPI::mpi.rank()) {
            // material from interior or copy
            auto neighbor = element.neighbors[side];
            initAssign(material.neighbor[side], materialsDB[neighbor]);
          } else {
            // material from ghost layer (computed locally)
            auto neighborRank = element.neighborRanks[side];
            auto neighborRankIdx = element.mpiIndices[side];
            auto materialGhostIdx = ghostIdxMap.at(neighborRank)[neighborRankIdx];
            initAssign(material.neighbor[side], materialsDBGhost[materialGhostIdx]);
          }
        } else {
          // otherwise, use the material from the own cell
          initAssign(material.neighbor[side], localMaterial);
        }
      }

      // if enabled, set up the plasticity as well
      if (seissolParams.model.plasticity) {
        auto& plasticity = plasticityArray[cell];
        const auto& localPlasticity = plasticityDB[meshId];

        initAssign(plasticity, seissol::model::PlasticityData(localPlasticity, &material.local));
      }
    }
    ltsToMesh += layer.getNumberOfCells();
  }
}

struct LtsInfo {
  unsigned* ltsMeshToFace = nullptr;
  MeshStructure* meshStructure = nullptr;
  TimeStepping timeStepping{};

  // IMPORTANT: DO NOT DEALLOCATE THE ABOVE POINTERS... THEY ARE PASSED ON AND REQUIRED DURING
  // RUNTIME
};

void initializeCellMatrices(LtsInfo& ltsInfo, seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();

  // \todo Move this to some common initialization place
  auto& meshReader = seissolInstance.meshReader();
  auto& memoryManager = seissolInstance.getMemoryManager();

  seissol::initializer::initializeCellLocalMatrices(meshReader,
                                                    memoryManager.getLtsTree(),
                                                    memoryManager.getLts(),
                                                    memoryManager.getLtsLut(),
                                                    ltsInfo.timeStepping,
                                                    seissolParams.model);

  if (seissolParams.drParameters.etaHack != 1.0) {
    logWarning(seissol::MPI::mpi.rank())
        << "The \"eta hack\" has been enabled to mitigate quasi-divergent solutions in the "
           "friction law. The results may not conform to the existing benchmarks.";
  }

  seissol::initializer::initializeDynamicRuptureMatrices(meshReader,
                                                         memoryManager.getLtsTree(),
                                                         memoryManager.getLts(),
                                                         memoryManager.getLtsLut(),
                                                         memoryManager.getDynamicRuptureTree(),
                                                         memoryManager.getDynamicRupture(),
                                                         ltsInfo.ltsMeshToFace,
                                                         *memoryManager.getGlobalDataOnHost(),
                                                         seissolParams.drParameters.etaHack);

  memoryManager.initFrictionData();

  seissol::initializer::initializeBoundaryMappings(meshReader,
                                                   memoryManager.getEasiBoundaryReader(),
                                                   memoryManager.getLtsTree(),
                                                   memoryManager.getLts(),
                                                   memoryManager.getLtsLut());

#ifdef ACL_DEVICE
  initializer::copyCellMatricesToDevice(memoryManager.getLtsTree(),
                                        memoryManager.getLts(),
                                        memoryManager.getDynamicRuptureTree(),
                                        memoryManager.getDynamicRupture(),
                                        memoryManager.getBoundaryTree(),
                                        memoryManager.getBoundary());

  memoryManager.recordExecutionPaths(seissolParams.model.plasticity);
#endif

  auto itmParameters = seissolInstance.getSeisSolParameters().model.itmParameters;

  if (itmParameters.itmEnabled) {
    auto& timeMirrorManagers = seissolInstance.getTimeMirrorManagers();
    const double scalingFactor = itmParameters.itmVelocityScalingFactor;
    const double startingTime = itmParameters.itmStartingTime;

    auto* mLtsTree = memoryManager.getLtsTree();
    auto* mLts = memoryManager.getLts();
    auto* mLtsLut = memoryManager.getLtsLut();
    const auto* mTimeStepping = seissolInstance.timeManager().getTimeStepping();

    initializeTimeMirrorManagers(scalingFactor,
                                 startingTime,
                                 &meshReader,
                                 mLtsTree,
                                 mLts,
                                 mLtsLut,
                                 timeMirrorManagers.first,
                                 timeMirrorManagers.second,
                                 seissolInstance,
                                 mTimeStepping);
  }
}

void initializeClusteredLts(LtsInfo& ltsInfo, seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();

  assert(seissolParams.timeStepping.lts.getRate() > 0);

  if (seissolParams.timeStepping.lts.getRate() == 1) {
    seissolInstance.getLtsLayout().deriveLayout(TimeClustering::Single, 1);
  } else {
    seissolInstance.getLtsLayout().deriveLayout(TimeClustering::MultiRate,
                                                seissolParams.timeStepping.lts.getRate());
  }

  seissolInstance.getLtsLayout().getMeshStructure(ltsInfo.meshStructure);
  seissolInstance.getLtsLayout().getCrossClusterTimeStepping(ltsInfo.timeStepping);

  seissolInstance.getMemoryManager().initializeFrictionLaw();

  unsigned* numberOfDRCopyFaces = nullptr;
  unsigned* numberOfDRInteriorFaces = nullptr;

  seissolInstance.getLtsLayout().getDynamicRuptureInformation(
      ltsInfo.ltsMeshToFace, numberOfDRCopyFaces, numberOfDRInteriorFaces);

  seissolInstance.getMemoryManager().fixateLtsTree(ltsInfo.timeStepping,
                                                   ltsInfo.meshStructure,
                                                   numberOfDRCopyFaces,
                                                   numberOfDRInteriorFaces,
                                                   seissolParams.model.plasticity);

  seissolInstance.getMemoryManager().setLtsToFace(ltsInfo.ltsMeshToFace);

  delete[] numberOfDRCopyFaces;
  delete[] numberOfDRInteriorFaces;

  const auto& ltsTree = seissolInstance.getMemoryManager().getLtsTree();
  const auto& lts = seissolInstance.getMemoryManager().getLts();

  unsigned* ltsToMesh = nullptr;
  unsigned numberOfMeshCells = 0;

  seissolInstance.getLtsLayout().getCellInformation(
      ltsTree->var(lts->cellInformation), ltsToMesh, numberOfMeshCells);

  // TODO(David): move all of this method to the MemoryManager
  seissolInstance.getMemoryManager().getLtsLutUnsafe().createLuts(
      ltsTree, ltsToMesh, numberOfMeshCells);

  delete[] ltsToMesh;

  seissol::initializer::time_stepping::deriveLtsSetups(ltsInfo.timeStepping.numberOfLocalClusters,
                                                       ltsInfo.meshStructure,
                                                       ltsTree->var(lts->cellInformation));
}

void initializeMemoryLayout(LtsInfo& ltsInfo, seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();

  seissolInstance.getMemoryManager().initializeMemoryLayout();

  seissolInstance.timeManager().addClusters(ltsInfo.timeStepping,
                                            ltsInfo.meshStructure,
                                            seissolInstance.getMemoryManager(),
                                            seissolParams.model.plasticity);

  // set tv for all time clusters (this needs to be done, after the time clusters start existing)
  if (seissolParams.model.plasticity) {
    seissolInstance.timeManager().setTv(seissolParams.model.tv);
  }

  seissolInstance.getMemoryManager().fixateBoundaryLtsTree();
}

} // namespace

void seissol::initializer::initprocedure::initModel(seissol::SeisSol& seissolInstance) {
  SCOREP_USER_REGION("init_model", SCOREP_USER_REGION_TYPE_FUNCTION);

  logInfo(seissol::MPI::mpi.rank()) << "Begin init model.";

  // Call the pre mesh initialization hook
  seissol::Modules::callHook<ModuleHook::PreModel>();

  seissol::Stopwatch watch;
  watch.start();

  LtsInfo ltsInfo;

  // these four methods need to be called in this order.

  // init LTS
  logInfo(seissol::MPI::mpi.rank()) << "Initialize LTS.";
  initializeClusteredLts(ltsInfo, seissolInstance);

  // init cell materials (needs LTS, to place the material in; this part was translated from
  // FORTRAN)
  logInfo(seissol::MPI::mpi.rank()) << "Initialize cell material parameters.";
  initializeCellMaterial(seissolInstance);

  // init memory layout (needs cell material values to initialize e.g. displacements correctly)
  logInfo(seissol::MPI::mpi.rank()) << "Initialize Memory layout.";
  initializeMemoryLayout(ltsInfo, seissolInstance);

  // init cell matrices
  logInfo(seissol::MPI::mpi.rank()) << "Initialize cell-local matrices.";
  initializeCellMatrices(ltsInfo, seissolInstance);

  watch.pause();
  watch.printTime("Model initialized in:");

  // Call the post mesh initialization hook
  seissol::Modules::callHook<ModuleHook::PostModel>();

  logInfo(seissol::MPI::mpi.rank()) << "End init model.";
}
