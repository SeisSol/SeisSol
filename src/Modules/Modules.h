/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 */

#ifndef MODULES_H
#define MODULES_H

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
  PreMPI,
  PostMPIInit,
  PreMesh,
  PostMesh,
  PreLtsInit,
  PostLtsInit,
  PreModel,
  PostModel,
  /**
   * Called when the simulation starts.
   *
   * @warning Only called when the simulation is not loaded from a checkpoint.
   */
  SimulationStart,
  /**
   * Global synchronization point during simulation
   *
   * Registering for this hook requires setting the update interval.
   */
  SynchronizationPoint,
  SimulationEnd,
  Shutdown,
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
  ModuleHook nextHook;

  private:
  Modules();

  /**
   * Do the real work for registering for a hook
   */
  void _registerHook(Module& module, ModuleHook hook, ModulePriority priority);

  /**
   * Do the real work for handling a hook
   */
  template <ModuleHook hook>
  void _callHook() {
    for (auto& [_, module] : hooks[static_cast<size_t>(ModuleHook::SynchronizationPoint)]) {
      call<hook>(module);
    }

    nextHook = static_cast<ModuleHook>(static_cast<int>(hook) + 1);
  }

  double _callSyncHook(double currentTime, double timeTolerance, bool forceSyncPoint);

  /**
   * Set the simulation start time.
   *
   * This is required to handle synchronization points correctly when the simulation starts
   * from a checkpoint.
   */
  void _setSimulationStartTime(double time);

  private:
  template <ModuleHook hook>
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

  template <ModuleHook hook>
  static void callHook() {
    instance()._callHook<hook>();
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

#endif // MODULES_H
