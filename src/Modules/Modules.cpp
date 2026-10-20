// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#include "Modules.h"

#include "Memory/HostData.h"
#include "Modules/Module.h"
#include "Modules/OutputWorker.h"

#include <algorithm>
#include <async/as/Pin.h>
#include <cassert>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol {

void Modules::_registerHook(Module& module, ModuleHook hook, ModulePriority priority) {
  assert(static_cast<int>(hook) < static_cast<int>(ModuleHook::MaxHooks));

  if (nextHook_ >= ModuleHook::MaxInitHooks) {
    logError() << "Trying to register for a hook after initialization phase";
  }
  if (hook < nextHook_) {
    logError() << "Trying to register for hook" << strHook(hook)
               << "but SeisSol was already processing"
               << strHook(static_cast<ModuleHook>(static_cast<int>(nextHook_) - 1));
  }

  hooks_[static_cast<size_t>(hook)].insert(std::pair<ModulePriority, Module*>(priority, &module));
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

Modules::~Modules() = default;

double Modules::_callSyncHook(double currentTime, double timeTolerance, bool forceSyncPoint) {
  auto& modules = hooks_[static_cast<size_t>(ModuleHook::SynchronizationPoint)];

  // With a thread for the outputs, the modules due here process the synchronization point there,
  // if all of them can, one after the other, as they would on the main thread; so they keep their
  // order (e.g. of their collectives) on all ranks.
  std::vector<Module*> due;
  for (auto& [_, module] : modules) {
    if (module->dueAt(currentTime, timeTolerance, forceSyncPoint)) {
      due.push_back(module);
    }
  }
  const bool asynchronous = worker_ != nullptr && !due.empty() &&
                            std::all_of(due.begin(), due.end(), [](const Module* module) {
                              return module->processesAsynchronously();
                            });

  double nextSyncTime = std::numeric_limits<double>::max();
  if (asynchronous) {
    for (auto* module : due) {
      module->captureSyncPoint(currentTime);
      module->markSyncPoint(currentTime);
    }
    worker_->enqueue([due, currentTime] {
      for (auto* module : due) {
        module->processSyncPoint(currentTime);
      }
    });
    for (auto& [_, module] : modules) {
      nextSyncTime = std::min(nextSyncTime, module->nextSyncPoint());
    }
  } else {
    if (worker_ != nullptr) {
      worker_->wait();
    }
    for (auto& [_, module] : modules) {
      nextSyncTime = std::min(
          nextSyncTime, module->potentialSyncPoint(currentTime, timeTolerance, forceSyncPoint));
    }
  }

  return nextSyncTime;
}

initializer::HostData
    Modules::_hostDataAt(double currentTime, double timeTolerance, bool forceSyncPoint) const {
  initializer::HostData data;
  for (const auto& [_, module] : hooks_[static_cast<size_t>(ModuleHook::SynchronizationPoint)]) {
    if (module->dueAt(currentTime, timeTolerance, forceSyncPoint)) {
      data.merge(module->hostData());
    }
  }
  return data;
}

void Modules::_callSimulationStartHook(std::optional<double> checkpointTime) {
  assert(static_cast<int>(nextHook_) <= static_cast<int>(ModuleHook::SynchronizationPoint));

  const auto startTime = checkpointTime.value_or(0);

  for (auto& [_, module] : hooks_[static_cast<size_t>(ModuleHook::SimulationStart)]) {
    module->simulationStart(checkpointTime);
  }

  // Modules that register only for synchronization points need a valid
  // initial sync schedule (e.g., ITM managers).
  for (auto& [_, module] : hooks_[static_cast<size_t>(ModuleHook::SynchronizationPoint)]) {
    module->setSimulationStartTime(startTime);
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

initializer::HostData
    Modules::hostDataAt(double currentTime, double timeTolerance, bool forceSyncPoint) {
  return instance()._hostDataAt(currentTime, timeTolerance, forceSyncPoint);
}

void Modules::startOutputWorker(const async::as::CpuMask& cpus, int ompThreads) {
  assert(instance().worker_ == nullptr);
  instance().worker_ = std::make_unique<OutputWorker>(cpus, ompThreads);
}

void Modules::waitForOutputs() {
  if (instance().worker_ != nullptr) {
    instance().worker_->wait();
  }
}

bool Modules::outputWorkerRunning() { return instance().worker_ != nullptr; }

void Modules::stopOutputWorker() {
  // the thread processes what is left before it ends
  instance().worker_.reset();
}

void Modules::callSimulationStartHook(std::optional<double> checkpointTime) {
  instance()._callSimulationStartHook(checkpointTime);
}

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
MODULES_CALL_INSTANCE(ModuleHook::SimulationEnd, simulationEnd)
MODULES_CALL_INSTANCE(ModuleHook::Shutdown, shutdown)

} // namespace seissol
