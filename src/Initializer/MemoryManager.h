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

#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
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
  //! global data of every configuration
  PerConfig<GlobalData> globalDataOnHost_;
  PerConfig<GlobalData> globalDataOnDevice_;

  //! Memory organization storage
  LTS::Storage ltsStorage_;
  LTS::Backmap backmap_;

  // the initial conditions of the cells of each configuration, by its id
  std::vector<std::vector<std::unique_ptr<physics::InitialField>>> iniConds_;

  DynamicRupture::Backmap drBackmap_;

  DynamicRupture::Storage drStorage_;
  std::unique_ptr<DynamicRupture> dynRup_ = nullptr;
  std::unique_ptr<dr::initializer::BaseDRInitializer> drInitializer_ = nullptr;
  dr::friction_law::FrictionSolverFactory frictionLaw_;
  dr::friction_law::FrictionSolverFactory frictionLawDevice_;
  std::unique_ptr<dr::output::OutputManager> faultOutputManager_ = nullptr;

  Boundary::Storage boundaryStorage_;

  SurfaceLTS::Storage surfaceStorage_;

  std::optional<ClusterLayout> layout_;

  public:
  /**
   * Default constructor
   **/
  explicit MemoryManager(seissol::SeisSol& instance);

  /**
   * Initialization function, which allocates memory for the global matrices and initializes them.
   **/
  void initialize();

  /**
   * Gets the global data of the configuration `Cfg` on both host and device.
   **/
  template <typename Cfg>
  CompoundGlobalData<Cfg> globalData() {
    CompoundGlobalData<Cfg> global{};
    global.onHost = &globalDataOnHost_.get<Cfg>();
    global.onDevice = nullptr;
    if constexpr (seissol::isDeviceOn()) {
      global.onDevice = &globalDataOnDevice_.get<Cfg>();
    }
    return global;
  }

  void setClusterLayout(const ClusterLayout& extLayout) { layout_.emplace(extLayout); }

  ClusterLayout& clusterLayout() { return layout_.value(); }

  LTS::Storage& ltsStorage() { return ltsStorage_; }

  LTS::Backmap& backmap() { return backmap_; }

  DynamicRupture::Storage& drStorage() { return drStorage_; }

  DynamicRupture::Backmap& drBackmap() { return drBackmap_; }

  DynamicRupture& drDescriptor() { return *dynRup_; }

  SurfaceLTS::Storage& surfaceStorage() { return surfaceStorage_; }

  Boundary::Storage& boundaryStorage() { return boundaryStorage_; }

  /// Sets the initial conditions of the cells of the configuration `config`.
  void setInitialConditions(ConfigId config,
                            std::vector<std::unique_ptr<physics::InitialField>>&& iniConds) {
    iniConds_.at(config) = std::move(iniConds);
  }

  /// The initial conditions of the cells of the configuration `config`; the reference stays valid.
  const std::vector<std::unique_ptr<physics::InitialField>>& initialConditions(ConfigId config) {
    return iniConds_.at(config);
  }

  /// The initial conditions of the cells of every configuration, by its id.
  const std::vector<std::vector<std::unique_ptr<physics::InitialField>>>& initialConditions() {
    return iniConds_;
  }

  const dr::friction_law::FrictionSolverFactory& frictionLaw() { return frictionLaw_; }
  const dr::friction_law::FrictionSolverFactory& frictionLawDevice() { return frictionLawDevice_; }
  seissol::dr::output::OutputManager* faultOutputManager() { return faultOutputManager_.get(); }

  void initializeFrictionLaw();
  void initFrictionData();
  void synchronizeTo(seissol::initializer::AllocationPlace place);
};
} // namespace initializer
} // namespace seissol

#endif // SEISSOL_SRC_INITIALIZER_MEMORYMANAGER_H_
