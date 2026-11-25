// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_SIMULATION_MEMORYMANAGER_H_
#define SEISSOL_SRC_SIMULATION_MEMORYMANAGER_H_

#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/InputAux.h"
#include "Core/Kernels/GlobalDataStore.h"
#include "Core/Memory/Descriptor/Boundary.h"
#include "Core/Memory/Descriptor/DynamicRupture.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Descriptor/Surface.h"
#include "Core/Memory/MemoryAllocator.h"
#include "Core/Memory/SimulationSettings.h"
#include "Core/Memory/Tree/Backmap.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Memory/Typedefs.h"
#include "Core/TimeStepping/Layout/ClusterLayout.h"
#include "Features/DynamicRupture/Setup/Factory.h"
#include "Simulation/Parameters/SeisSolParameters.h"
#include "WavePropagation/Physics/InitialField.h"

#include <memory>
#include <mpi.h>
#include <optional>
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
  GlobalDataStore globalData_;

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

  std::optional<SimulationSettings> simulationSettings_;

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
    return globalData_.get<Cfg>();
  }

  /// The global data of all configurations.
  GlobalDataStore& globalDataStore() { return globalData_; }

  void setClusterLayout(const ClusterLayout& extLayout) { layout_.emplace(extLayout); }

  ClusterLayout& clusterLayout() { return layout_.value(); }

  void setSimulationSettings(const SimulationSettings& settings) {
    simulationSettings_.emplace(settings);
  }

  /// The settings of the simulations, as the cells are laid out for them.
  [[nodiscard]] const SimulationSettings& simulationSettings() const {
    return simulationSettings_.value();
  }

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

  /**
   * Copies the given variables to the host, and waits for that; all variables of all storages for
   * everything.
   */
};
} // namespace initializer
} // namespace seissol

#endif // SEISSOL_SRC_SIMULATION_MEMORYMANAGER_H_
