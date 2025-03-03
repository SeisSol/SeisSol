// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_MODULES_MODULES_H_
#define SEISSOL_SRC_MODULES_MODULES_H_

#include "Module.h"

#include "utils/logger.h"
#include <array>
#include <limits>
#include <map>

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
  PreMPI = 0,
  PostMPIInit = 1,
  PreMesh = 2,
  PostMesh = 3,
  PreLtsInit = 4,
  PostLtsInit = 5,
  PreModel = 6,
  PostModel = 7,
  /**
   * Called when the simulation starts.
   *
   * @warning Only called when the simulation is not loaded from a checkpoint.
   */
  SimulationStart = 8,
  /**
   * Global synchronization point during simulation
   *
   * Registering for this hook requires setting the update interval.
   */
  SynchronizationPoint = 9,
  SimulationEnd = 10,
  Shutdown = 11,
  FirstHook = PreMPI,
  MaxInitHooks = SimulationStart + 1,
  MaxHooks = Shutdown + 1
};

/**
 * Manages SeisSol modules
 */
class Modules {
  private:
  std::array<std::multimap<ModulePriority, Module*>, static_cast<size_t>(ModuleHook::MaxHooks)>
      hooks;

  /** The hook that should be called next */
  ModuleHook nextHook{ModuleHook::FirstHook};

  Modules();

  /**
   * Do the real work for registering for a hook
   */
  void _registerHook(Module& module, ModuleHook hook, ModulePriority priority);

  /**
   * Do the real work for handling a hook
   */
  template <ModuleHook Hook>
  void _callHook() {
    for (auto& [_, module] : hooks[static_cast<size_t>(ModuleHook::SynchronizationPoint)]) {
      call<Hook>(module);
    }

    nextHook = static_cast<ModuleHook>(static_cast<int>(Hook) + 1);
  }

  double _callSyncHook(double currentTime, double timeTolerance, bool forceSyncPoint);

  /**
   * Set the simulation start time.
   *
   * This is required to handle synchronization points correctly when the simulation starts
   * from a checkpoint.
   */
  void _setSimulationStartTime(double time);

  template <ModuleHook Hook>
  static void call(Module* module);

  static const char* strHook(ModuleHook hook);

  /**
   * The only instance of this class
   *
   * We need to use a static function to prevent the "static initialization order fiasco".
   */
  static Modules& instance();

  // Public interface
  public:
  static void registerHook(Module& module,
                           ModuleHook hook,
                           ModulePriority priority = ModulePriority::Default);

  template <ModuleHook Hook>
  static void callHook() {
    instance()._callHook<Hook>();
  }

  /**
   * @param currentTime The current simulation time.
   * @param timeTolerance The time tolerance for time comparison
   * @return The next synchronization point
   *
   * @todo The time tolerance is global constant, maybe not necessary to pass it here
   */
  static double callSyncHook(double currentTime, double timeTolerance, bool forceSyncPoint = false);

  /**
   * Set the simulation start time
   */
  static void setSimulationStartTime(double time);
};

template <>
inline void seissol::Modules::_callHook<ModuleHook::SynchronizationPoint>() {
  logError() << "Synchronization point hooks have to be called with \"callSyncHook\"";
}

} // namespace seissol

#endif // SEISSOL_SRC_MODULES_MODULES_H_
