// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#include "Modules/Module.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <utility>
#include <utils/logger.h>

#include "Modules.h"

namespace seissol {

void Modules::_registerHook(Module& module, ModuleHook hook, ModulePriority priority) {
  assert(static_cast<int>(hook) < static_cast<int>(ModuleHook::MaxHooks));

  if (nextHook >= ModuleHook::MaxInitHooks) {
    logError() << "Trying to register for a hook after initialization phase";
  }
  if (hook < nextHook) {
    logError() << "Trying to register for hook" << strHook(hook)
               << "but SeisSol was already processing"
               << strHook(static_cast<ModuleHook>(static_cast<int>(nextHook) - 1));
  }

  hooks[static_cast<size_t>(hook)].insert(std::pair<ModulePriority, Module*>(priority, &module));
}

const char* Modules::strHook(ModuleHook hook) {
  switch (hook) {
  case ModuleHook::PreMPI:
    return "PRE_MPI";
  case ModuleHook::PostMPIInit:
    return "POST_MPI_INIT";
  case ModuleHook::PreMesh:
    return "PRE_MESH";
  case ModuleHook::PostMesh:
    return "POST_MESH";
  case ModuleHook::PreLtsInit:
    return "PRE_LTSINIT";
  case ModuleHook::PostLtsInit:
    return "POST_LTSINIT";
  case ModuleHook::PreModel:
    return "PRE_MODEL";
  case ModuleHook::PostModel:
    return "POST_MODEL";
  case ModuleHook::SimulationStart:
    return "SIMULATION_START";
  case ModuleHook::SynchronizationPoint:
    return "SYNCHRONIZATION_POINT";
  case ModuleHook::SimulationEnd:
    return "SIMULATION_END";
  case ModuleHook::Shutdown:
    return "SHUTDOWN";
  default:
    return "unknown hook";
  }
}

Modules::Modules() = default;

double Modules::_callSyncHook(double currentTime, double timeTolerance, bool forceSyncPoint) {
  double nextSyncTime = std::numeric_limits<double>::max();

  for (auto& [_, module] : hooks[static_cast<size_t>(ModuleHook::SynchronizationPoint)]) {
    nextSyncTime = std::min(nextSyncTime,
                            module->potentialSyncPoint(currentTime, timeTolerance, forceSyncPoint));
  }

  return nextSyncTime;
}

void Modules::_setSimulationStartTime(double time) {
  assert(static_cast<int>(nextHook) <= static_cast<int>(ModuleHook::SynchronizationPoint));

  // Set the simulation time in all modules that are called at synchronization points
  for (auto& [_, module] : hooks[static_cast<size_t>(ModuleHook::SynchronizationPoint)]) {
    module->setSimulationStartTime(time);
  }
}

Modules& Modules::instance() {
  static Modules instance;
  return instance;
}

void Modules::registerHook(Module& module, ModuleHook hook, ModulePriority priority) {
  instance()._registerHook(module, hook, priority);
}

double Modules::callSyncHook(double currentTime, double timeTolerance, bool forceSyncPoint) {
  return instance()._callSyncHook(currentTime, timeTolerance, forceSyncPoint);
}

void Modules::setSimulationStartTime(double time) { instance()._setSimulationStartTime(time); }

// Create all template instances for call
#define MODULES_CALL_INSTANCE(enum, func)                                                          \
  template <>                                                                                      \
  void Modules::call<enum>(Module * module) {                                                      \
    module->func();                                                                                \
  }

MODULES_CALL_INSTANCE(ModuleHook::PreMPI, preMPI)
MODULES_CALL_INSTANCE(ModuleHook::PostMPIInit, postMPIInit)
MODULES_CALL_INSTANCE(ModuleHook::PreMesh, preMesh)
MODULES_CALL_INSTANCE(ModuleHook::PostMesh, postMesh)
MODULES_CALL_INSTANCE(ModuleHook::PreLtsInit, preLtsInit)
MODULES_CALL_INSTANCE(ModuleHook::PostLtsInit, postLtsInit)
MODULES_CALL_INSTANCE(ModuleHook::PreModel, preModel)
MODULES_CALL_INSTANCE(ModuleHook::PostModel, postModel)
MODULES_CALL_INSTANCE(ModuleHook::SimulationStart, simulationStart)
MODULES_CALL_INSTANCE(ModuleHook::SimulationEnd, simulationEnd)
MODULES_CALL_INSTANCE(ModuleHook::Shutdown, shutdown)

} // namespace seissol
