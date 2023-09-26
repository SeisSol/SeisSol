
#include <Initializer/ConfigFile.hpp>
#include <vector>
#include "Initializer/ParameterDB.h"
#include "Initializer/InputParameters.hpp"
#include "Initializer/CellLocalMatrices.h"
#include "Initializer/LTS.h"
#include "Initializer/tree/LTSTree.hpp"
#include "Physics/Attenuation.hpp"
#include "Equations/datastructures.hpp"
#include "Initializer/tree/Lut.hpp"
#include "Initializer/typedefs.hpp"
#include "Initializer/ParameterMaterialDB.hpp"

#include "SeisSol.h"
#include "Init.hpp"
#include "InitModel.hpp"

#include "Parallel/MPI.h"

#include "Initializer/Layout/Memory.hpp"

#include <cmath>
#include <type_traits>

using namespace seissol::initializer;

using Plasticity = seissol::model::Plasticity;

/*
Assigns the given value to the target object, initializing the memory in the process.

NOTE: std::copy (or the likes) do not work here, since they do not initialize the _vptr for virtual
function calls (rather, they leave it undefined), since they do merely assign `value` to `target`.
*/

template <typename T>
static void initAssign(T& target, const T& value) {
  if constexpr (std::is_trivially_copyable_v<T>) {
    // if the object is trivially copyable, we may just memcpy it (it's safe to do that in this
    // case).
    std::memcpy(&target, &value, sizeof(T));
  } else {
    // otherwise, call the class/struct initializer.
    // problem: we may have an array here... So we unwrap it.
    if constexpr (std::is_array_v<T>) {
      // unwrap array, dimension by dimension...
      // example: T[N][M] yields SubT=T[M]
      using SubT = std::remove_extent_t<T>;
      auto subExtent = std::extent_v<T>;

      // for now, init element-wise... (TODO(David): we could look for something faster here, in
      // case it should ever matter)
      for (size_t i = 0; i < subExtent; ++i) {
        initAssign<SubT>(target[i], value[i]);
      }
    } else {
      // now call new here.
      new (&target) T(value);
    }
  }
  // (these two methods cannot be combined, unless we have some way for C-style arrays, i.e. S[N]
  // for <typename S, size_t N>, to use a copy constructor as well)
}

template <typename T>
static void synchronize(const seissol::initializers::Variable<T>& handle) {
  auto& memoryManager = seissol::SeisSol::main.getMemoryManager();
  const auto& meshToLts = memoryManager.getLtsLut()->getMeshToLtsLut(handle.mask);
  unsigned* duplicatedMeshIds = memoryManager.getLtsLut()->getDuplicatedMeshIds(handle.mask);
  const unsigned numberOfDuplicatedMeshIds =
      memoryManager.getLtsLut()->getNumberOfDuplicatedMeshIds(handle.mask);
  T* var = ltsview.tree.var(handle);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (unsigned dupMeshId = 0; dupMeshId < numberOfDuplicatedMeshIds; ++dupMeshId) {
    const unsigned meshId = duplicatedMeshIds[dupMeshId];
    const T* ref = &var[meshToLts[0][meshId]];
    for (unsigned dup = 1; dup < seissol::initializers::Lut::MaxDuplicates &&
                           meshToLts[dup][meshId] != std::numeric_limits<unsigned>::max();
         ++dup) {

      // copy data on a byte-wise level (we need to initialize memory here as well)
      initAssign(var[meshToLts[dup][meshId]], *ref);
    }
  }
}

/*void initializeCellMaterial() {
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

  const configs = seissol::initializer::readConfigFile(seissolParams.model.configFileName);

  // material retrieval for interior/copy and ghost layers
  auto [materialsDB, plasticityDB] = seissol::model::queryMaterial(
      configs, seissol::initializers::CellToVertexArray::fromMeshReader(meshReader), true);
  auto [materialsDBGhost, plasticityDBGhost] = seissol::model::queryMaterial(
      configs,
      seissol::initializers::CellToVertexArray::fromVectors(ghostVertices, ghostGroups),
      true);

  logDebug()
      << "Setting cell materials in the LTS tree."; // TODO(David): describe plasticity as well
  const auto& elements = meshReader.getElements();

  unsigned* ltsToMesh =
      memoryManager.getLtsLut()->getLtsToMeshLut(seissol::initializers::LayerMask(Ghost));

  memoryManager.getClusterForest().visit([&](auto& ltsview) {
    auto* materialDataGlobal = ltsview.tree.var(ltsview.lts.materialData);

    for (seissol::initializers::LTSTree::leaf_iterator it =
             ltsview.tree.beginLeaf(seissol::initializers::LayerMask(Ghost));
         it != ltsview.tree.endLeaf();
         ++it) {
      auto* cellInformation = it->var(ltsview.lts.cellInformation);
      auto* materialLayer = it->var(ltsview.lts.material);
      auto* materialDataLayer = it->var(ltsview.lts.materialData);
      auto* plasticityArray =
          seissolParams.model.plasticity ? it->var(ltsview.lts.plasticity) : nullptr;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (std::size_t cell = 0; cell < it->getNumberOfCells(); ++cell) {
        // this loop does three things at the same time:

        // set the materials for the cell volume and its faces
        auto meshId = ltsToMesh[cell];
        auto& material = materialLayer[cell];
        auto& materialData = materialDataLayer[cell];
        const auto& element = elements[meshId];
        const auto& localCellInformation = cellInformation[cell];

        // write material data
        initAssign(materialData, materialsDB[meshId]);
        for (std::size_t side = 0; side < 4; ++side) {
          // set ghost rank materials
          if (element.neighborRanks[side] != seissol::MPI::mpi.rank()) {
            // material from ghost layer (computed locally)
            auto neighborId = localCellInformation.faceNeighborIds[side];

            auto neighborRank = element.neighborRanks[side];
            auto neighborRankIdx = element.mpiIndices[side];
            auto materialGhostIdx = ghostIdxMap.at(neighborRank)[neighborRankIdx];
            initAssign(materialDataGlobal[neighborId], materialsDBGhost[materialGhostIdx]);
          }
          // all other elements are set up locally in some other iteration
        }

        // set material pointers
        material.local = &materialData;
        for (std::size_t side = 0; side < 4; ++side) {
          auto neighborId = localCellInformation.faceNeighborIds[side];
          if (isInternalFaceType(localCellInformation.faceTypes[side])) {
            // use neighbor material for internal/interior faces
            material.neighbor[side] = &materialDataGlobal[neighborId];
          } else {
            // otherwise, use the material from the own cell
            material.neighbor[side] = &materialData;
          }
        }

        // if enabled, set up the plasticity as well
        if constexpr (seissolParams.model.plasticity) {
          plasticityArray[cell] =
              seissol::model::PlasticityData<RealT>(plasticityDB[meshId], material.local);
        }
      }
      ltsToMesh += it->getNumberOfCells();
    }
  });
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
                                                     memoryManager.getClusterForest(),
                                                     memoryManager.getLtsLut(),
                                                     ltsInfo.timeStepping);

  seissol::initializers::initializeDynamicRuptureMatrices(meshReader,
                                                          memoryManager.getClusterForest(),
                                                          memoryManager.getLtsLut(),
                                                          memoryManager.getDynamicRuptureForest(),
                                                          ltsInfo.ltsMeshToFace,
                                                          *memoryManager.getGlobalDataOnHost(),
                                                          ltsInfo.timeStepping);

  memoryManager.initFrictionData();

  seissol::initializers::initializeBoundaryMappings(meshReader,
                                                    memoryManager.getEasiBoundaryReader(),
                                                    memoryManager.getClusterForest(),
                                                    memoryManager.getLtsLut());

#ifdef ACL_DEVICE
  initializers::copyCellMatricesToDevice(memoryManager.getClusterForest(),
                                         memoryManager.getDynamicRuptureForest(),
                                         memoryManager.getBoundaryForest());

  memoryManager.recordExecutionPaths(seissolParams.model.plasticity);
#endif

  // synchronize data
  synchronize(ltsview.lts.dofs);
  if (kernels::size<tensor::Qane>() > 0) {
    synchronize(ltsview.lts.dofsAne);
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

  // set tv for all time clusters
  if (seissolParams.model.plasticity) {
    seissol::SeisSol::main.timeManager().setTv(seissolParams.model.tv);
  }

  seissol::SeisSol::main.getMemoryManager().fixateBoundaryLtsTree();
}*/

void seissol::initializer::initprocedure::initModel() {
  SCOREP_USER_REGION("init_memory", SCOREP_USER_REGION_TYPE_FUNCTION);

  logInfo(seissol::MPI::mpi.rank()) << "Begin init memory.";

  // Call the pre mesh initialization hook
  seissol::Modules::callHook<seissol::PRE_MODEL>();

  seissol::Stopwatch watch;
  watch.start();

  /*LtsInfo ltsInfo;

  // these four methods need to be called in this order.

  // init LTS (needed in all subsequent steps)
  logInfo(seissol::MPI::mpi.rank()) << "Initialize LTS.";
  initializeClusteredLts(ltsInfo);

  // init cell materials (needs LTS, to place the material in)
  logInfo(seissol::MPI::mpi.rank()) << "Initialize cell material parameters.";
  initializeCellMaterial();

  // init memory layout (needs cell material values to initialize e.g. displacements correctly)
  logInfo(seissol::MPI::mpi.rank()) << "Initialize memory layout.";
  initializeMemoryLayout(ltsInfo);

  // init cell matrices (needs the data from all previous steps)
  logInfo(seissol::MPI::mpi.rank()) << "Initialize cell-local matrices.";
  initializeCellMatrices(ltsInfo);*/

  auto memoryContainer =
      setupMemory(seissol::SeisSol::main.meshReader(), seissol::SeisSol::main.timeManager());

  // seissol::SeisSol::main.getMemoryManager().setMemoryContainer(std::move(memoryContainer));

  watch.pause();
  watch.printTime("Memory initialized in:");

  // Call the post mesh initialization hook
  seissol::Modules::callHook<seissol::POST_MODEL>();

  logInfo(seissol::MPI::mpi.rank()) << "End init memory.";
}
