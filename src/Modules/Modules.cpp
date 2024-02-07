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

#include <cassert>

#include "Modules.h"

namespace seissol {

void Modules::_registerHook(Module& module, ModuleHook hook, ModulePriority priority) {
  assert(hook < MAX_HOOKS);

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

Modules::Modules() : nextHook(ModuleHook::FirstHook) {}

double Modules::_callSyncHook(double currentTime, double timeTolerance, bool forceSyncPoint) {
  double nextSyncTime = std::numeric_limits<double>::max();

  for (auto& [_, module] : hooks[static_cast<size_t>(ModuleHook::SynchronizationPoint)]) {
    nextSyncTime = std::min(nextSyncTime,
                            module->potentialSyncPoint(currentTime, timeTolerance, forceSyncPoint));
  }

  return nextSyncTime;
}

void Modules::_setSimulationStartTime(double time) {
  assert(nextHook <= static_cast<size_t>(ModuleHook::SynchronizationPoint));

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

void Modules::setSimulationStartTime(double time) {
  return instance()._setSimulationStartTime(time);
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
MODULES_CALL_INSTANCE(ModuleHook::SimulationStart, simulationStart)
MODULES_CALL_INSTANCE(ModuleHook::SimulationEnd, simulationEnd)
MODULES_CALL_INSTANCE(ModuleHook::Shutdown, shutdown)

} // namespace seissol
