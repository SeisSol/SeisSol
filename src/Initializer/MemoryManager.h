// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_INITIALIZER_MEMORYMANAGER_H_
#define SEISSOL_SRC_INITIALIZER_MEMORYMANAGER_H_

#include "Initializer/Parameters/SeisSolParameters.h"
#include "Memory/Tree/Layer.h"
#include <Config.h>
#include <Initializer/TimeStepping/ClusterLayout.h>
#include <Memory/MemoryContainer.h>
#include <Memory/Tree/Backmap.h>
#include <mpi.h>

#include <utils/logger.h>

#include "Initializer/Typedefs.h"
#include "Memory/MemoryAllocator.h"

#include "Initializer/InputAux.h"
#include "Initializer/ParameterDB.h"
#include "Memory/Descriptor/Boundary.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/LTSTree.h"

#include "Physics/InitialField.h"

#include <memory>
#include <utility>
#include <vector>

#include "DynamicRupture/Factory.h"
#include <yaml-cpp/yaml.h>

namespace seissol {
class SeisSol;
namespace initializer {

/**
 * Memory manager of SeisSol.
 **/
class MemoryManager {
  private: // explicit private for unit tests
  seissol::SeisSol& seissolInstance;

  std::optional<memory::MemoryContainer> container;

  //! memory allocator
  seissol::memory::ManagedAllocator m_memoryAllocator;

  //! LTS mesh structure
  struct MeshStructure* m_meshStructure;

  unsigned int* ltsToFace;

  /*
   * Cross-cluster
   */
  //! global data
  GlobalData m_globalDataOnHost;
  GlobalData m_globalDataOnDevice;

  //! Memory organization tree
  LTSTree m_ltsTree;
  LTS m_lts;

  std::vector<std::unique_ptr<physics::InitialField>> m_iniConds;

  LTSTree m_dynRupTree;
  std::unique_ptr<DynamicRupture> m_dynRup = nullptr;
  std::unique_ptr<dr::initializer::BaseDRInitializer> m_DRInitializer = nullptr;
  std::unique_ptr<dr::friction_law::FrictionSolver> m_FrictionLaw = nullptr;
  std::unique_ptr<dr::friction_law::FrictionSolver> m_FrictionLawDevice = nullptr;
  std::unique_ptr<dr::output::OutputManager> m_faultOutputManager = nullptr;
  std::shared_ptr<seissol::initializer::parameters::SeisSolParameters> m_seissolParams = nullptr;

  LTSTree m_boundaryTree;
  Boundary m_boundary;

  EasiBoundary m_easiBoundary;

  public:
  /**
   * Constructor
   **/
  MemoryManager(seissol::SeisSol& instance) : seissolInstance(instance) {};

  /**
   * Destructor, memory is freed by managed allocator
   **/
  ~MemoryManager() = default;

  /**
   * Initialization function, which allocates memory for the global matrices and initializes them.
   **/
  void initialize();

  void setupMemoryContainer(std::size_t maxCluster, const std::vector<ConfigVariant>& configs) {
    container.emplace(maxCluster, configs);
    container.value().globalDataStorage.onHost = &m_globalDataOnHost;
    container.value().globalDataStorage.onDevice = &m_globalDataOnDevice;
  }

  memory::MemoryContainer& memoryContainer() { return container.value(); }

  /**
   * Sets the number of cells in each leaf of the lts tree, fixates the variables, and allocates
   *memory. Afterwards the tree cannot be changed anymore.
   *
   * @param i_meshStructrue mesh structure.
   **/
  void fixateLtsTree(struct ClusterLayout& clusterLayout,
                     struct MeshStructure* meshStructure,
                     unsigned* numberOfDRCopyFaces,
                     unsigned* numberOfDRInteriorFaces,
                     bool usePlasticity);

  void fixateBoundaryLtsTree();
  /**
   * Set up the internal structure.
   **/
  void initializeMemoryLayout();

  /**
   * Gets the global data on both host and device.
   **/
  CompoundGlobalData getGlobalData() { return memoryContainer().globalDataStorage; }

  /**
   * Gets the memory layout of a time cluster.
   *
   * @param i_cluster local id of the time cluster.
   * @param o_meshStructure mesh structure.
   * @param o_globalData global data.
   * @param o_globalDataCopies several copies of global data
   **/
  std::pair<MeshStructure*, CompoundGlobalData> getMemoryLayout(unsigned int cluster);

  LTSTree* getLtsTree() { return &container.value().volume; }

  LTS* getLts() { return &container.value().wpdesc; }

  LTSTree* getDynamicRuptureTree() { return &container.value().dynrup; }

  DynamicRupture* getDynamicRupture() { return container.value().drdesc.get(); }

  LTSTree* getBoundaryTree() { return &container.value().boundary; }

  Boundary* getBoundary() { return &container.value().bnddesc; }

  StorageBackmap<4, LTSTree>* wpBackmap() { return &container.value().clusterBackmap; }

  StorageBackmap<1, LTSTree>* drBackmap() { return &container.value().dynrupBackmap; }

  void setInitialConditions(std::vector<std::unique_ptr<physics::InitialField>>&& iniConds) {
    m_iniConds = std::move(iniConds);
  }

  const std::vector<std::unique_ptr<physics::InitialField>>& getInitialConditions() {
    return m_iniConds;
  }

  void setLtsToFace(unsigned int* ptr) { ltsToFace = ptr; }

  unsigned int* ltsToFaceMap() const { return ltsToFace; }

  void initializeEasiBoundaryReader(const char* fileName);

  EasiBoundary* getEasiBoundaryReader() { return &m_easiBoundary; }

  dr::friction_law::FrictionSolver* getFrictionLaw() { return m_FrictionLaw.get(); }
  dr::friction_law::FrictionSolver* getFrictionLawDevice() { return m_FrictionLawDevice.get(); }
  dr::initializer::BaseDRInitializer* getDRInitializer() { return m_DRInitializer.get(); }
  seissol::dr::output::OutputManager* getFaultOutputManager() { return m_faultOutputManager.get(); }
  seissol::initializer::parameters::DRParameters* getDRParameters() {
    return &(m_seissolParams->drParameters);
  }

  seissol::initializer::parameters::LtsParameters* getLtsParameters() {
    return &(m_seissolParams->timeStepping.lts);
  };

  void setInputParams(std::shared_ptr<seissol::initializer::parameters::SeisSolParameters> params) {
    m_seissolParams = std::move(params);
  }

  std::string getOutputPrefix() const { return m_seissolParams->output.prefix; }

  bool isLoopStatisticsNetcdfOutputOn() const {
    return m_seissolParams->output.loopStatisticsNetcdfOutput;
  }

#ifdef ACL_DEVICE
  void recordExecutionPaths(bool usePlasticity);

  /**
   * Derives sizes of scratch memory required during computations of Wave Propagation solver
   **/
  static void deriveRequiredScratchpadMemoryForWp(bool plasticity, LTSTree& ltsTree, LTS& lts);

  /**
   * Derives sizes of scratch memory required during computations of Dynamic Rupture solver
   **/
  static void deriveRequiredScratchpadMemoryForDr(LTSTree& ltsTree, DynamicRupture& dynRup);
#endif

  void initializeFrictionLaw();
  void initFaultOutputManager(const std::string& backupTimeStamp);
  void initFrictionData();
  void synchronizeTo(seissol::initializer::AllocationPlace place);
};

bool isAcousticSideOfElasticAcousticInterface(const CellMaterialData& material, unsigned int face);
bool isElasticSideOfElasticAcousticInterface(const CellMaterialData& material, unsigned int face);
bool isAtElasticAcousticInterface(const CellMaterialData& material, unsigned int face);

bool requiresDisplacement(CellLocalInformation cellLocalInformation,
                          const CellMaterialData& material,
                          unsigned int face);
bool requiresNodalFlux(FaceType f);
} // namespace initializer
} // namespace seissol

#endif // SEISSOL_SRC_INITIALIZER_MEMORYMANAGER_H_
