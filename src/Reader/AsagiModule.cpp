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
 * Velocity field reader Fortran interface
 */

#ifdef USE_ASAGI

#include "AsagiModule.h"
#include "utils/env.h"
#include "Parallel/Helper.hpp"
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

namespace seissol::asagi {

AsagiModule::AsagiModule() : m_mpiMode(getMPIMode()), m_totalThreads(getTotalThreads()) {
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

AsagiMPIMode AsagiModule::getMPIMode() {
#ifdef USE_MPI
  std::string mpiModeName = utils::Env::get(EnvMPIMode, "WINDOWS");
  if (mpiModeName == "WINDOWS")
    return AsagiMPIMode::Windows;
  if (mpiModeName == "COMM_THREAD")
    return AsagiMPIMode::CommThread;
  if (mpiModeName == "OFF")
    return AsagiMPIMode::Off;

  return AsagiMPIMode::Unknown;
#else  // USE_MPI
  return AsagiMPIMode::Off;
#endif // USE_MPI
}

int AsagiModule::getTotalThreads() {
  int totalThreads = 1;

#ifdef _OPENMP
  totalThreads = omp_get_max_threads();
  if (seissol::useCommThread(seissol::MPI::mpi)) {
    totalThreads++;
  }
#endif // _OPENMP

  return totalThreads;
}

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
    std::string mpiModeName = utils::Env::get(EnvMPIMode, "");
    logError() << "Unknown ASAGI MPI mode:" << mpiModeName;
  } else {
    const int rank = MPI::mpi.rank();
    logWarning(rank) << "Running with only one OMP thread."
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

AsagiModule& AsagiModule::getInstance() {
  static AsagiModule instance;
  return instance;
}

AsagiMPIMode AsagiModule::mpiMode() { return getInstance().m_mpiMode; }

int AsagiModule::totalThreads() { return getInstance().m_totalThreads; }

} // namespace seissol::asagi

#endif
