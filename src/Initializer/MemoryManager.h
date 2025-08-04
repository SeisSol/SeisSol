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
#include <Memory/Descriptor/Surface.h>
#include <Memory/GlobalData.h>
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

  //! memory allocator
  seissol::memory::ManagedAllocator m_memoryAllocator;

  /*
   * Cross-cluster
   */
  //! global data
  GlobalData global;

  //! Memory organization storage
  LTS::Storage ltsStorage;
  LTS::Backmap backmap;

  std::vector<std::unique_ptr<physics::InitialField>> m_iniConds;

  DynamicRupture::Backmap drBackmap;

  DynamicRupture::Storage drStorage;
  std::unique_ptr<DynamicRupture> m_dynRup = nullptr;
  std::unique_ptr<dr::initializer::BaseDRInitializer> m_DRInitializer = nullptr;
  dr::friction_law::FrictionSolverFactory m_FrictionLaw = nullptr;
  dr::friction_law::FrictionSolverFactory m_FrictionLawDevice = nullptr;
  std::unique_ptr<dr::output::OutputManager> m_faultOutputManager = nullptr;

  Boundary::Storage m_boundaryTree;

  SurfaceLTS::Storage surfaceStorage;

  EasiBoundaryT m_easiBoundary;

  std::optional<ClusterLayout> layout;

  public:
  /**
   * Constructor
   **/
  MemoryManager(seissol::SeisSol& instance) : seissolInstance(instance) {}

  /**
   * Destructor, memory is freed by managed allocator
   **/
  ~MemoryManager() = default;

  /**
   * Initialization function, which allocates memory for the global matrices and initializes them.
   **/
  void initialize();

  /**
   * Sets the number of cells in each leaf of the lts storage, fixates the variables, and allocates
   *memory. Afterwards the storage cannot be changed anymore.
   *
   * @param i_meshStructrue mesh structure.
   **/
  void fixateLtsStorage();

  void fixateBoundaryStorage();
  /**
   * Set up the internal structure.
   **/
  void initializeMemoryLayout();

  /**
   * Gets the global data on both host and device.
   **/
  GlobalData& getGlobalData() { return global; }

  void setClusterLayout(const ClusterLayout& extLayout) { layout.emplace(extLayout); }

  ClusterLayout& clusterLayout() { return layout.value(); }

  LTS::Storage& getLtsStorage() { return ltsStorage; }

  LTS::Backmap& getBackmap() { return backmap; }

  DynamicRupture::Storage& getDRStorage() { return drStorage; }

  DynamicRupture::Backmap& getDRBackmap() { return drBackmap; }

  DynamicRupture& getDynamicRupture() { return *m_dynRup; }

  SurfaceLTS::Storage& getSurfaceStorage() { return surfaceStorage; }

  void setInitialConditions(std::vector<std::unique_ptr<physics::InitialField>>&& iniConds) {
    m_iniConds = std::move(iniConds);
  }

  const std::vector<std::unique_ptr<physics::InitialField>>& getInitialConditions() {
    return m_iniConds;
  }

  void initializeEasiBoundaryReader(const char* fileName);

  EasiBoundary* getEasiBoundaryReader() { return &m_easiBoundary; }

  dr::friction_law::FrictionSolverFactory& getFrictionLaw() { return m_FrictionLaw; }
  dr::friction_law::FrictionSolverFactory& getFrictionLawDevice() { return m_FrictionLawDevice; }
  seissol::dr::output::OutputManager* getFaultOutputManager() { return m_faultOutputManager.get(); }

#ifdef ACL_DEVICE
  void recordExecutionPaths(bool usePlasticity);

  /**
   * Derives sizes of scratch memory required during computations of Wave Propagation solver
   **/
  static void deriveRequiredScratchpadMemoryForWp(bool plasticity, LTS::Storage& ltsStorage);

  /**
   * Derives sizes of scratch memory required during computations of Dynamic Rupture solver
   **/
  static void deriveRequiredScratchpadMemoryForDr(DynamicRupture::Storage& drStorage);
#endif

  void initializeFrictionLaw();
  void initFaultOutputManager(const std::string& backupTimeStamp);
  void initFrictionData();
  void synchronizeTo(seissol::initializer::AllocationPlace place);
};

bool isAcousticSideOfElasticAcousticInterface(CellMaterialData& material, std::size_t face);
bool isElasticSideOfElasticAcousticInterface(CellMaterialData& material, std::size_t face);
bool isAtElasticAcousticInterface(CellMaterialData& material, std::size_t face);

bool requiresDisplacement(CellLocalInformation cellLocalInformation,
                          CellMaterialData& material,
                          std::size_t face);
bool requiresNodalFlux(FaceType f);
} // namespace initializer
} // namespace seissol

#endif // SEISSOL_SRC_INITIALIZER_MEMORYMANAGER_H_
