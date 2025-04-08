// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifdef USE_ASAGI

#include "AsagiModule.h"
#include "Parallel/Helper.h"
#include "utils/env.h"
#include <Modules/Modules.h>
#include <Parallel/MPI.h>
#include <asagi.h>
#include <memory>
#include <string>
#include <utils/logger.h>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

namespace seissol::asagi {

AsagiModule::AsagiModule(utils::Env& env)
    : m_env(env), m_mpiMode(getMPIMode(env)), m_totalThreads(getTotalThreads(env)) {
  // Register for the pre MPI hook
  Modules::registerHook(*this, ModuleHook::PreMPI);

  // Emit a warning/error later
  // TODO use a general logger that can buffer log messages and emit them later
  if (m_mpiMode == AsagiMPIMode::Unknown) {
    Modules::registerHook(*this, ModuleHook::PostMPIInit);
  } else if (m_mpiMode == AsagiMPIMode::CommThread && m_totalThreads == 1) {
    m_mpiMode = AsagiMPIMode::Windows;

    Modules::registerHook(*this, ModuleHook::PostMPIInit);
  }
}

AsagiMPIMode AsagiModule::getMPIMode(utils::Env& env) {
#ifdef USE_MPI
  const std::string mpiModeName = env.get(EnvMpiMode, "WINDOWS");
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
  int totalThreads = 1;

#ifdef _OPENMP
  totalThreads = omp_get_max_threads();
  if (seissol::useCommThread(seissol::MPI::mpi, env)) {
    totalThreads++;
  }
#endif // _OPENMP

  return totalThreads;
}

utils::Env& AsagiModule::getEnv() { return m_env; }

void AsagiModule::preMPI() {
  // Communication threads required
  if (m_mpiMode == AsagiMPIMode::CommThread) {
    // Comm threads has to be started before model initialization
    Modules::registerHook(*this, ModuleHook::PreModel, ModulePriority::Highest);
    // Comm threads has to be stoped after model initialization
    Modules::registerHook(*this, ModuleHook::PostModel, ModulePriority::Lowest);
  }
}

void AsagiModule::postMPIInit() {
  if (m_mpiMode == AsagiMPIMode::Unknown) {
    const std::string mpiModeName = m_env.get(EnvMpiMode, "");
    logError() << "Unknown ASAGI MPI mode:" << mpiModeName;
  } else {
    logWarning() << "Running with only one OMP thread."
                 << "Using MPI window communication instead of threads.";
  }
}

void AsagiModule::preModel() {
#ifdef USE_MPI
  // TODO check if ASAGI is required for model setup
  ::asagi::Grid::startCommThread();
#endif // USE_MPI
}

void AsagiModule::postModel() {
#ifdef USE_MPI
  // TODO check if ASAGI is required for model setup
  ::asagi::Grid::stopCommThread();
#endif // USE_MPI
}

void AsagiModule::initInstance(utils::Env& env) {
  AsagiModule::instance = std::make_shared<AsagiModule>(env);
}

AsagiModule& AsagiModule::getInstance() { return *AsagiModule::instance; }

AsagiMPIMode AsagiModule::mpiMode() { return getInstance().m_mpiMode; }

int AsagiModule::totalThreads() { return getInstance().m_totalThreads; }

std::shared_ptr<AsagiModule> AsagiModule::instance{nullptr};

} // namespace seissol::asagi

#endif
