
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
#include "Initializer/tree/LTSSync.hpp"
#include "Initializer/typedefs.hpp"

#include "SeisSol.h"
#include "Init.hpp"
#include "InitModel.hpp"

#include "Parallel/MPI.h"

#include <cmath>
#include <type_traits>

using namespace seissol::initializer;

namespace {

using Material_t = seissol::model::Material_t;
using Plasticity = seissol::model::Plasticity;

template <typename T>
static std::vector<T> queryDB(seissol::initializers::QueryGenerator* queryGen,
                              const std::string& fileName,
                              size_t size) {
  std::vector<T> vectorDB(size);
  seissol::initializers::MaterialParameterDB<T> parameterDB;
  parameterDB.setMaterialVector(&vectorDB);
  parameterDB.evaluateModel(fileName, queryGen);
  return vectorDB;
}

void initializeCellMaterial() {
  const auto& seissolParams = seissol::SeisSol::main.getSeisSolParameters();
  const auto& meshReader = seissol::SeisSol::main.meshReader();
  initializers::MemoryManager& memoryManager = seissol::SeisSol::main.getMemoryManager();

  // unpack ghost layer (merely a re-ordering operation, since the CellToVertexArray right now
  // requires an vector there)
  std::vector<std::array<std::array<double, 3>, 4>> ghostVertices;
  std::vector<int> ghostGroups;
  std::unordered_map<int, std::vector<unsigned>> ghostIdxMap;
  for (const auto& neighbor : meshReader.getGhostlayerMetadata()) {
    ghostIdxMap[neighbor.first].reserve(neighbor.second.size());
    for (const auto& metadata : neighbor.second) {
      ghostIdxMap[neighbor.first].push_back(ghostVertices.size());
      std::array<std::array<double, 3>, 4> vertices;
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
  auto getBestQueryGenerator = [&](const seissol::initializers::CellToVertexArray& ctvArray) {
    return seissol::initializers::getBestQueryGenerator(
        seissol::initializer::parameters::isModelAnelastic(),
        seissolParams.model.plasticity,
        seissol::initializer::parameters::isModelAnisotropic(),
        seissol::initializer::parameters::isModelPoroelastic(),
        seissolParams.model.useCellHomogenizedMaterial,
        ctvArray);
  };

  // material retrieval for copy+interior layers
  seissol::initializers::QueryGenerator* queryGen =
      getBestQueryGenerator(seissol::initializers::CellToVertexArray::fromMeshReader(meshReader));
  auto materialsDB = queryDB<Material_t>(
      queryGen, seissolParams.model.materialFileName, meshReader.getElements().size());

  // plasticity (if needed)
  std::vector<Plasticity> plasticityDB;
  if (seissolParams.model.plasticity) {
    // plasticity information is only needed on all interior+copy cells.
    plasticityDB = queryDB<Plasticity>(
        queryGen, seissolParams.model.materialFileName, meshReader.getElements().size());
  }

  // material retrieval for ghost layers
  seissol::initializers::QueryGenerator* queryGenGhost = getBestQueryGenerator(
      seissol::initializers::CellToVertexArray::fromVectors(ghostVertices, ghostGroups));
  auto materialsDBGhost = queryDB<Material_t>(
      queryGenGhost, seissolParams.model.materialFileName, ghostVertices.size());

#if defined(USE_VISCOELASTIC) || defined(USE_VISCOELASTIC2)
  // we need to compute all model parameters before we can use them...
  // TODO(David): integrate this with the Viscoelastic material class or the ParameterDB directly?
  logDebug() << "Initializing attenuation.";
#ifdef OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (size_t i = 0; i < materialsDB.size(); ++i) {
    auto& cellMat = materialsDB[i];
    seissol::physics::fitAttenuation(
        cellMat, seissolParams.model.freqCentral, seissolParams.model.freqRatio);
  }
#ifdef OPENMP
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
  unsigned* ltsToMesh =
      memoryManager.getLtsLut()->getLtsToMeshLut(memoryManager.getLts()->material.mask);

  for (seissol::initializers::LTSTree::leaf_iterator it =
           memoryManager.getLtsTree()->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != memoryManager.getLtsTree()->endLeaf();
       ++it) {
    auto* cellInformation = it->var(memoryManager.getLts()->cellInformation);
    auto* materialArray = it->var(memoryManager.getLts()->material);
    auto* plasticityArray =
        seissolParams.model.plasticity ? it->var(memoryManager.getLts()->plasticity) : nullptr;

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
      // TODO(David): move to material initalization maybe? Or an initializer for the PlasticityData
      // struct?
      if (seissolParams.model.plasticity) {
        auto& plasticity = plasticityArray[cell];
        const auto& localPlasticity = plasticityDB[meshId];

        plasticity.initialLoading[0] = localPlasticity.s_xx;
        plasticity.initialLoading[1] = localPlasticity.s_yy;
        plasticity.initialLoading[2] = localPlasticity.s_zz;
        plasticity.initialLoading[3] = localPlasticity.s_xy;
        plasticity.initialLoading[4] = localPlasticity.s_yz;
        plasticity.initialLoading[5] = localPlasticity.s_xz;

        const double angularFriction = std::atan(localPlasticity.bulkFriction);

        plasticity.cohesionTimesCosAngularFriction =
            localPlasticity.plastCo * std::cos(angularFriction);
        plasticity.sinAngularFriction = std::sin(angularFriction);
#ifndef USE_ANISOTROPIC
        plasticity.mufactor = 1.0 / (2.0 * material.local.mu);
#else
        plasticity.mufactor =
            3.0 / (2.0 * (material.local.c44 + material.local.c55 + material.local.c66));
#endif
      }
    }
    ltsToMesh += it->getNumberOfCells();
  }
}

struct LtsInfo {
  unsigned* ltsMeshToFace = nullptr;
  MeshStructure* meshStructure = nullptr;
  TimeStepping timeStepping;

  // IMPORTANT: DO NOT DEALLOCATE THE ABOVE POINTERS... THEY ARE PASSED ON AND REQUIRED DURING
  // RUNTIME
};

static void initializeCellMatrices(LtsInfo& ltsInfo) {
  const auto& seissolParams = seissol::SeisSol::main.getSeisSolParameters();

  // \todo Move this to some common initialization place
  auto& meshReader = seissol::SeisSol::main.meshReader();
  auto& memoryManager = seissol::SeisSol::main.getMemoryManager();

  seissol::initializers::initializeCellLocalMatrices(meshReader,
                                                     memoryManager.getLtsTree(),
                                                     memoryManager.getLts(),
                                                     memoryManager.getLtsLut(),
                                                     ltsInfo.timeStepping);

  seissol::initializers::initializeDynamicRuptureMatrices(meshReader,
                                                          memoryManager.getLtsTree(),
                                                          memoryManager.getLts(),
                                                          memoryManager.getLtsLut(),
                                                          memoryManager.getDynamicRuptureTree(),
                                                          memoryManager.getDynamicRupture(),
                                                          ltsInfo.ltsMeshToFace,
                                                          *memoryManager.getGlobalDataOnHost(),
                                                          ltsInfo.timeStepping);

  memoryManager.initFrictionData();

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

  memoryManager.recordExecutionPaths(seissolParams.model.plasticity);
#endif

  auto itmParameters = seissol::SeisSol::main.getSeisSolParameters().itmParameters;

  bool ITMToggle = itmParameters.ITMToggle;

  if (ITMToggle) {
    auto& timeMirrorManagers = seissol::SeisSol::main.getTimeMirrorManagers();
    double scalingFactor = itmParameters.ITMVelocityScalingFactor;
    double startingTime = itmParameters.ITMStartingTime;

    auto m_ltsTree = memoryManager.getLtsTree();
    auto m_lts = memoryManager.getLts();
    auto m_ltsLut = memoryManager.getLtsLut();
    auto m_timeStepping = seissol::SeisSol::main.timeManager().getTimeStepping();

    initializeTimeMirrorManagers(scalingFactor,
                                 startingTime,
                                 &meshReader,
                                 m_ltsTree,
                                 m_lts,
                                 m_ltsLut,
                                 timeMirrorManagers.first,
                                 timeMirrorManagers.second,
                                 m_timeStepping);
  }
}

static void initializeClusteredLts(LtsInfo& ltsInfo) {
  const auto& seissolParams = seissol::SeisSol::main.getSeisSolParameters();

  assert(seissolParams.timeStepping.lts.rate > 0);

  if (seissolParams.timeStepping.lts.rate == 1) {
    seissol::SeisSol::main.getLtsLayout().deriveLayout(single, 1);
  } else {
    seissol::SeisSol::main.getLtsLayout().deriveLayout(multiRate,
                                                       seissolParams.timeStepping.lts.rate);
  }

  seissol::SeisSol::main.getLtsLayout().getMeshStructure(ltsInfo.meshStructure);
  seissol::SeisSol::main.getLtsLayout().getCrossClusterTimeStepping(ltsInfo.timeStepping);

  seissol::SeisSol::main.getMemoryManager().initializeFrictionLaw();

  unsigned* numberOfDRCopyFaces;
  unsigned* numberOfDRInteriorFaces;

  seissol::SeisSol::main.getLtsLayout().getDynamicRuptureInformation(
      ltsInfo.ltsMeshToFace, numberOfDRCopyFaces, numberOfDRInteriorFaces);

  seissol::SeisSol::main.getMemoryManager().fixateLtsTree(ltsInfo.timeStepping,
                                                          ltsInfo.meshStructure,
                                                          numberOfDRCopyFaces,
                                                          numberOfDRInteriorFaces,
                                                          seissolParams.model.plasticity);

  delete[] numberOfDRCopyFaces;
  delete[] numberOfDRInteriorFaces;

  const auto& ltsTree = seissol::SeisSol::main.getMemoryManager().getLtsTree();
  const auto& lts = seissol::SeisSol::main.getMemoryManager().getLts();

  unsigned* ltsToMesh;
  unsigned numberOfMeshCells;

  seissol::SeisSol::main.getLtsLayout().getCellInformation(
      ltsTree->var(lts->cellInformation), ltsToMesh, numberOfMeshCells);

  // TODO(David): move all of this method to the MemoryManager
  seissol::SeisSol::main.getMemoryManager().getLtsLutUnsafe().createLuts(
      ltsTree, ltsToMesh, numberOfMeshCells);

  delete[] ltsToMesh;

  seissol::initializers::time_stepping::deriveLtsSetups(ltsInfo.timeStepping.numberOfLocalClusters,
                                                        ltsInfo.meshStructure,
                                                        ltsTree->var(lts->cellInformation));
}

static void initializeMemoryLayout(LtsInfo& ltsInfo) {
  const auto& seissolParams = seissol::SeisSol::main.getSeisSolParameters();

  seissol::SeisSol::main.getMemoryManager().initializeMemoryLayout();

  seissol::SeisSol::main.timeManager().addClusters(ltsInfo.timeStepping,
                                                   ltsInfo.meshStructure,
                                                   seissol::SeisSol::main.getMemoryManager(),
                                                   seissolParams.model.plasticity);

  // set tv for all time clusters (this needs to be done, after the time clusters start existing)
  if (seissolParams.model.plasticity) {
    seissol::SeisSol::main.timeManager().setTv(seissolParams.model.tv);
  }

  seissol::SeisSol::main.getMemoryManager().fixateBoundaryLtsTree();
}

} // namespace

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

  // init cell materials (needs LTS, to place the material in; this part was translated from
  // FORTRAN)
  logInfo(seissol::MPI::mpi.rank()) << "Initialize cell material parameters.";
  initializeCellMaterial();

  // init memory layout (needs cell material values to initialize e.g. displacements correctly)
  logInfo(seissol::MPI::mpi.rank()) << "Initialize Memory layout.";
  initializeMemoryLayout(ltsInfo);

  // init cell matrices
  logInfo(seissol::MPI::mpi.rank()) << "Initialize cell-local matrices.";
  initializeCellMatrices(ltsInfo);

  watch.pause();
  watch.printTime("Model initialized in:");

  // Call the post mesh initialization hook
  seissol::Modules::callHook<seissol::POST_MODEL>();

  logInfo(seissol::MPI::mpi.rank()) << "End init model.";
}
