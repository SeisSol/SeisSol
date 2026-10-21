// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_MODULES_MODULES_H_
#define SEISSOL_SRC_MODULES_MODULES_H_

#include "Module.h"

#include <array>
#include <limits>
#include <map>
#include <optional>
#include <utils/logger.h>

namespace seissol {

enum class ModulePriority : int {
  Max = std::numeric_limits<int>::min(),
  Highest = -10000,
  Higher = -1000,
  High = -100,
  Default = 0,
  Low = 100,
  Lower = 1000,
  Lowest = 10000,
  Min = std::numeric_limits<int>::max()
};

/**
 * Possible hooks for modules
 *
 * To add a new hook, add an additional enum entry here, add the corresponding
 * function to {@link Module}, add a specialization of {@link Modules::call(Module*)}
 * and update the function {@link Modules::strHook}.
 *
 * @warning The order of the hooks has to be the same they are called in SeisSol.
 */
enum class ModuleHook : int {
  NullHook = -1,
  PostMPIInit = 0,
  PreMesh = 1,
  PostMesh = 2,
  PreLtsInit = 3,
  PostLtsInit = 4,
  PreModel = 5,
  PostModel = 6,
  /**
   * Called when the simulation starts.
   *
   * @warning Only called when the simulation is not loaded from a checkpoint.
   */
  SimulationStart = 7,
  /**
   * Global synchronization point during simulation
   *
   * Registering for this hook requires setting the update interval.
   */
  SynchronizationPoint = 8,
  SimulationEnd = 9,
  Shutdown = 10,
  FirstHook = PostMPIInit,
  MaxInitHooks = SimulationStart + 1,
  MaxHooks = Shutdown + 1
};

/**
 * Manages the modules of a single SeisSol instance.
 *
 * The registry holds non-owning pointers to the modules that register with it, so it has to
 * outlive all of them.
 */
class Modules {
  public:
  Modules();

  /**
   * Register a module for a hook.
   *
   * Only permitted before the initialization phase has passed the hook in question.
   */
  void registerHook(Module& module,
                    ModuleHook hook,
                    ModulePriority priority = ModulePriority::Default);

  template <ModuleHook Hook>
  void callHook() {
    for (auto& [_, module] : hooks_[static_cast<size_t>(Hook)]) {
      call<Hook>(module);
    }

    nextHook_ = static_cast<ModuleHook>(static_cast<int>(Hook) + 1);
  }

  /**
   * @param currentTime The current simulation time.
   * @param timeTolerance The time tolerance for time comparison
   * @return The next synchronization point
   *
   * @todo The time tolerance is global constant, maybe not necessary to pass it here
   */
  double callSyncHook(double currentTime, double timeTolerance, bool forceSyncPoint = false);

  /**
   * Set the simulation start time.
   *
   * This is required to handle synchronization points correctly when the simulation starts
   * from a checkpoint.
   */
  void callSimulationStartHook(std::optional<double> checkpointTime);

  private:
  std::array<std::multimap<ModulePriority, Module*>, static_cast<size_t>(ModuleHook::MaxHooks)>
      hooks_;

  /** The hook that should be called next */
  ModuleHook nextHook_{ModuleHook::FirstHook};

  template <ModuleHook Hook>
  static void call(Module* module);

  static const char* strHook(ModuleHook hook);
};

template <>
inline void seissol::Modules::callHook<ModuleHook::SimulationStart>() {
  logError() << "Simulation start hooks have to be called with \"callSimulationStartHook\"";
}

template <>
inline void seissol::Modules::callHook<ModuleHook::SynchronizationPoint>() {
  logError() << "Synchronization point hooks have to be called with \"callSyncHook\"";
}

} // namespace seissol

#endif // SEISSOL_SRC_MODULES_MODULES_H_
