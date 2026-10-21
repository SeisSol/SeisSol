// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifdef USE_ASAGI

#include "AsagiModule.h"

#include "Modules/Modules.h"
#include "Parallel/Helper.h"
#include "Parallel/MPI.h"
#include "Parallel/OpenMP.h"
#include "Parallel/Pin.h"

#include <asagi.h>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <utils/env.h>
#include <utils/logger.h>

namespace seissol::asagi {

AsagiModule::AsagiModule(utils::Env& env)
    : env_(env), mpiMode_(getMPIMode(env)), totalThreads_(getTotalThreads(env)) {
  // Register for the pre MPI hook
  Modules::registerHook(*this, ModuleHook::PreMPI);

  // Emit a warning/error later
  // TODO use a general logger that can buffer log messages and emit them later
  if (mpiMode_ == AsagiMPIMode::Unknown) {
    Modules::registerHook(*this, ModuleHook::PostMPIInit);
  } else if (mpiMode_ == AsagiMPIMode::CommThread && totalThreads_ == 1) {
    mpiMode_ = AsagiMPIMode::Windows;

    Modules::registerHook(*this, ModuleHook::PostMPIInit);
  }
}

AsagiMPIMode AsagiModule::getMPIMode(utils::Env& env) {
  // USE_MPI kept on purpose
#ifdef USE_MPI
  const std::string mpiModeName = env.get(EnvMpiMode, "OFF");
  if (mpiModeName == "WINDOWS") {
    return AsagiMPIMode::Windows;
  }
  if (mpiModeName == "COMM_THREAD") {
    return AsagiMPIMode::CommThread;
  }
  if (mpiModeName == "OFF") {
    return AsagiMPIMode::Off;
  }

  return AsagiMPIMode::Unknown;
#else  // USE_MPI
  return AsagiMPIMode::Off;
#endif // USE_MPI
}

int AsagiModule::getTotalThreads(utils::Env& env) {
  int totalThreads = OpenMP::threadCount();
  if (seissol::useCommThread(seissol::Mpi::mpi, env)) {
    totalThreads++;
  }

  return totalThreads;
}

utils::Env& AsagiModule::getEnv() { return env_; }

void AsagiModule::setPinning(const parallel::Pinning& pinning) { pinning_ = &pinning; }

void AsagiModule::preMPI() {
  // Communication threads required
  if (mpiMode_ == AsagiMPIMode::CommThread) {
    // Comm threads has to be started before model initialization
    Modules::registerHook(*this, ModuleHook::PreMesh, ModulePriority::Highest);
    // Comm threads has to be stoped after model initialization
    Modules::registerHook(*this, ModuleHook::PostModel, ModulePriority::Lowest);
  }
}

void AsagiModule::postMPIInit() {
  if (mpiMode_ == AsagiMPIMode::Unknown) {
    const std::string mpiModeName = env_.get(EnvMpiMode, "");
    logError() << "Unknown ASAGI MPI mode:" << mpiModeName;
  } else {
    logWarning() << "Running with only one OMP thread."
                 << "Using MPI window communication instead of threads.";
  }
}

void AsagiModule::preMesh() {
  if (mpiMode_ == AsagiMPIMode::CommThread) {
    // ASAGI leaves the thread unpinned if the CPU id is not valid; a negative id counts from the
    // last CPU, so -1 is a valid one
    int cpu = std::numeric_limits<int>::max();
#ifndef __APPLE__
    if (pinning_ == nullptr) {
      logError() << "The ASAGI communication thread needs the pinning of the SeisSol instance.";
    } else {
      // (collective on the ranks of the node)
      const auto freeCpu = parallel::Pinning::firstCpu(pinning_->getFreeCPUsMask());
      if (freeCpu.has_value()) {
        cpu = freeCpu.value();
        logInfo() << "ASAGI communication thread affinity: CPU" << cpu;
      } else {
        logWarning() << "No free CPU left for the ASAGI communication thread; it is not pinned.";
      }
    }
#endif
    ::asagi::Grid::startCommThread(cpu);
  }
}

void AsagiModule::postModel() {
  if (mpiMode_ == AsagiMPIMode::CommThread) {
    ::asagi::Grid::stopCommThread();
  }
}

void AsagiModule::initInstance(utils::Env& env) {
  AsagiModule::instance = std::make_shared<AsagiModule>(env);
}

AsagiModule& AsagiModule::getInstance() { return *AsagiModule::instance; }

AsagiMPIMode AsagiModule::mpiMode() { return getInstance().mpiMode_; }

int AsagiModule::totalThreads() { return getInstance().totalThreads_; }

std::shared_ptr<AsagiModule> AsagiModule::instance{nullptr};

} // namespace seissol::asagi

#endif
