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

#include "DynamicRupture/Factory.h"
#include "Initializer/InputAux.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/Boundary.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Descriptor/Surface.h"
#include "Memory/MemoryAllocator.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/Layer.h"
#include "Physics/InitialField.h"

#include <memory>
#include <mpi.h>
#include <utility>
#include <utils/logger.h>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol {
class SeisSol;
namespace initializer {

/**
 * Memory manager of SeisSol.
 **/
class MemoryManager {
  private: // explicit private for unit tests
  seissol::SeisSol& seissolInstance_;

  //! memory allocator
  seissol::memory::ManagedAllocator memoryAllocator_;

  /*
   * Cross-cluster
   */
  //! global data
  GlobalData globalDataOnHost_;
  GlobalData globalDataOnDevice_;

  //! Memory organization storage
  LTS::Storage ltsStorage_;
  LTS::Backmap backmap_;

  std::vector<std::unique_ptr<physics::InitialField>> iniConds_;

  DynamicRupture::Backmap drBackmap_;

  DynamicRupture::Storage drStorage_;
  std::unique_ptr<DynamicRupture> dynRup_ = nullptr;
  std::unique_ptr<dr::initializer::BaseDRInitializer> DRInitializer_ = nullptr;
  std::unique_ptr<dr::friction_law::FrictionSolver> FrictionLaw_ = nullptr;
  std::unique_ptr<dr::friction_law::FrictionSolver> FrictionLawDevice_ = nullptr;
  std::unique_ptr<dr::output::OutputManager> faultOutputManager_ = nullptr;

  Boundary::Storage boundaryTree_;

  SurfaceLTS::Storage surfaceStorage_;

  EasiBoundary easiBoundary_;

  std::optional<ClusterLayout> layout_;

  public:
  /**
   * Constructor
   **/
  explicit MemoryManager(seissol::SeisSol& instance) : seissolInstance_(instance) {}

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
  CompoundGlobalData getGlobalData() {
    CompoundGlobalData global{};
    global.onHost = &globalDataOnHost_;
    global.onDevice = nullptr;
    if constexpr (seissol::isDeviceOn()) {
      global.onDevice = &globalDataOnDevice_;
    }
    return global;
  }

  void setClusterLayout(const ClusterLayout& extLayout) { layout_.emplace(extLayout); }

  ClusterLayout& clusterLayout() { return layout_.value(); }

  LTS::Storage& getLtsStorage() { return ltsStorage_; }

  LTS::Backmap& getBackmap() { return backmap_; }

  DynamicRupture::Storage& getDRStorage() { return drStorage_; }

  DynamicRupture::Backmap& getDRBackmap() { return drBackmap_; }

  DynamicRupture& getDynamicRupture() { return *dynRup_; }

  SurfaceLTS::Storage& getSurfaceStorage() { return surfaceStorage_; }

  void setInitialConditions(std::vector<std::unique_ptr<physics::InitialField>>&& iniConds) {
    iniConds_ = std::move(iniConds);
  }

  const std::vector<std::unique_ptr<physics::InitialField>>& getInitialConditions() {
    return iniConds_;
  }

  void initializeEasiBoundaryReader(const char* fileName);

  EasiBoundary* getEasiBoundaryReader() { return &easiBoundary_; }

  dr::friction_law::FrictionSolver* getFrictionLaw() { return FrictionLaw_.get(); }
  dr::friction_law::FrictionSolver* getFrictionLawDevice() { return FrictionLawDevice_.get(); }
  seissol::dr::output::OutputManager* getFaultOutputManager() { return faultOutputManager_.get(); }

  void recordExecutionPaths(bool usePlasticity);

  /**
   * Derives sizes of scratch memory required during computations of Wave Propagation solver
   **/
  static void deriveRequiredScratchpadMemoryForWp(bool plasticity, LTS::Storage& ltsStorage);

  /**
   * Derives sizes of scratch memory required during computations of Dynamic Rupture solver
   **/
  static void deriveRequiredScratchpadMemoryForDr(DynamicRupture::Storage& drStorage);

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
