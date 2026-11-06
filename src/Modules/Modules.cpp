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
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol {

void Modules::registerHook(Module& module, ModuleHook hook, ModulePriority priority) {
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

std::vector<Module*>
    Modules::dueAt(double currentTime, double timeTolerance, bool forceSyncPoint) const {
  std::vector<Module*> due;
  for (const auto& [_, module] : hooks_[static_cast<size_t>(ModuleHook::SynchronizationPoint)]) {
    if (module->dueAt(currentTime, timeTolerance, forceSyncPoint)) {
      due.push_back(module);
    }
  }
  return due;
}

bool Modules::asynchronous(const std::vector<Module*>& due) const {
  return worker_ != nullptr && !due.empty() &&
         std::all_of(due.begin(), due.end(), [](const Module* module) {
           return module->processesAsynchronously();
         });
}

double Modules::callSyncHook(double currentTime,
                             double timeTolerance,
                             bool forceSyncPoint,
                             const std::function<void()>& awaitHostData) {
  auto& modules = hooks_[static_cast<size_t>(ModuleHook::SynchronizationPoint)];

  // With a thread for the outputs, the modules due here process the synchronization point there,
  // if all of them can, one after the other, as they would on the main thread; so they keep their
  // order (e.g. of their collectives) on all ranks.
  const auto due = dueAt(currentTime, timeTolerance, forceSyncPoint);

  double nextSyncTime = std::numeric_limits<double>::max();
  if (asynchronous(due)) {
    for (auto* module : due) {
      module->captureSyncPoint(currentTime);
      module->markSyncPoint(currentTime);
    }
    worker_->enqueue([due, currentTime, awaitHostData] {
      if (awaitHostData) {
        awaitHostData();
      }
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
    if (awaitHostData) {
      awaitHostData();
    }
    for (auto& [_, module] : modules) {
      nextSyncTime = std::min(
          nextSyncTime, module->potentialSyncPoint(currentTime, timeTolerance, forceSyncPoint));
    }
  }

  return nextSyncTime;
}

bool Modules::processesAsynchronouslyAt(double currentTime,
                                        double timeTolerance,
                                        bool forceSyncPoint) const {
  return asynchronous(dueAt(currentTime, timeTolerance, forceSyncPoint));
}

initializer::HostData
    Modules::hostDataAt(double currentTime, double timeTolerance, bool forceSyncPoint) const {
  initializer::HostData data;
  for (const auto& [_, module] : hooks_[static_cast<size_t>(ModuleHook::SynchronizationPoint)]) {
    if (module->dueAt(currentTime, timeTolerance, forceSyncPoint)) {
      data.merge(module->hostData());
    }
  }
  return data;
}

void Modules::startOutputWorker(const async::as::CpuMask& cpus, int ompThreads) {
  assert(worker_ == nullptr);
  worker_ = std::make_unique<OutputWorker>(cpus, ompThreads);
}

void Modules::waitForOutputs() {
  if (worker_ != nullptr) {
    worker_->wait();
  }
}

bool Modules::outputWorkerRunning() const { return worker_ != nullptr; }

void Modules::enqueueOnOutputThread(std::function<void()> task) {
  assert(worker_ != nullptr);
  worker_->enqueue(std::move(task));
}

void Modules::stopOutputWorker() {
  // the thread processes what is left before it ends
  worker_.reset();
}

void Modules::callSimulationStartHook(std::optional<double> checkpointTime) {
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

// Create all template instances for call
#define MODULES_CALL_INSTANCE(enum, func)                                                          \
  template <>                                                                                      \
  void Modules::call<enum>(Module * module) {                                                      \
    module->func();                                                                                \
  }

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
