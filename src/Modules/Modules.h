// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_MODULES_MODULES_H_
#define SEISSOL_SRC_MODULES_MODULES_H_

#include "Memory/HostData.h"
#include "Module.h"

#include <array>
#include <async/as/Pin.h>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <utils/logger.h>

namespace seissol {

class OutputWorker;

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
      hooks_;

  /** The hook that should be called next */
  ModuleHook nextHook_{ModuleHook::FirstHook};

  /** The thread the outputs may process their synchronization points on */
  std::unique_ptr<OutputWorker> worker_;

  Modules();
  ~Modules();

  /**
   * Do the real work for registering for a hook
   */
  void _registerHook(Module& module, ModuleHook hook, ModulePriority priority);

  /**
   * Do the real work for handling a hook
   */
  template <ModuleHook Hook>
  void _callHook() {
    for (auto& [_, module] : hooks_[static_cast<size_t>(Hook)]) {
      call<Hook>(module);
    }

    nextHook_ = static_cast<ModuleHook>(static_cast<int>(Hook) + 1);
  }

  double _callSyncHook(double currentTime, double timeTolerance, bool forceSyncPoint);

  [[nodiscard]] initializer::HostData
      _hostDataAt(double currentTime, double timeTolerance, bool forceSyncPoint) const;

  /**
   * Set the simulation start time.
   *
   * This is required to handle synchronization points correctly when the simulation starts
   * from a checkpoint.
   */
  void _callSimulationStartHook(std::optional<double> checkpointTime);

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
   * The variables of the storages that the modules due at a synchronization point read on the
   * host there; for the same arguments as callSyncHook(), before it.
   */
  [[nodiscard]] static initializer::HostData
      hostDataAt(double currentTime, double timeTolerance, bool forceSyncPoint = false);

  /**
   * Starts a thread on `cpus`, on which the modules due at a synchronization point process it
   * (see Module::processSyncPoint()) while the time stepping goes on, if all of them can; with up
   * to `ompThreads` threads in its OpenMP regions.
   */
  static void startOutputWorker(const async::as::CpuMask& cpus, int ompThreads);

  /** Waits until the synchronization points so far are processed. */
  static void waitForOutputs();

  /** Whether the outputs have a thread of their own. */
  [[nodiscard]] static bool outputWorkerRunning();

  /** Waits for the outputs, and ends the thread they are processed on. */
  static void stopOutputWorker();

  /**
   * Set the simulation start time
   */
  static void callSimulationStartHook(std::optional<double> checkpointTime);
};

template <>
inline void seissol::Modules::_callHook<ModuleHook::SimulationStart>() {
  logError() << "Simulation start hooks have to be called with \"callSimulationStartHook\"";
}

template <>
inline void seissol::Modules::_callHook<ModuleHook::SynchronizationPoint>() {
  logError() << "Synchronization point hooks have to be called with \"callSyncHook\"";
}

} // namespace seissol

#endif // SEISSOL_SRC_MODULES_MODULES_H_
