// SPDX-FileCopyrightText: 2014 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#include "SeisSol.h"

#include "Modules/Modules.h"
#include "Monitoring/Unit.h"
#include "Parallel/Helper.h"
#include "Parallel/MPI.h"
#include "Parallel/OpenMP.h"
#include "Parallel/Pin.h"

#include <cstddef>
#include <memory>
#include <optional>
#include <sys/resource.h>
#include <utils/logger.h>

namespace seissol {

bool SeisSol::init() {
  const auto rank = seissol::Mpi::mpi.rank();

  if (rank == 0) {
    const auto& hostNames = seissol::Mpi::mpi.getHostNames();
    logInfo() << "Running on (rank=0):" << hostNames.front();
  }

  logInfo() << "Using MPI with #ranks:" << seissol::Mpi::mpi.size();
  logInfo() << "Node-wide (shared memory) MPI with #ranks/node:"
            << seissol::Mpi::mpi.sharedMemMpiSize();
  seissol::Mpi::mpi.printAcceleratorDeviceInfo();
  // TODO (Ravil, David): switch to reading MPI options from the parameter-file.
  seissol::Mpi::mpi.setDataTransferModeFromEnv();

  printPersistentMpiInfo(env_);
#ifdef ACL_DEVICE
  printUSMInfo(env_);
  printMPIUSMInfo(env_);
#endif
  pinning_.checkEnvVariables();
  if (OpenMP::enabled()) {
    logInfo() << "Using OpenMP with #threads/rank:" << seissol::OpenMP::threadCount();
  } else {
    logInfo() << "OpenMP disabled. Using only a single thread.";
  }
  if (!parallel::Pinning::areAllCpusOnline()) {
    logInfo() << "Some CPUs are offline. Only online CPUs are considered.";
    logInfo() << "Online Mask            (this node)   :"
              << parallel::Pinning::maskToString(pinning_.getOnlineMask());
  }
  logInfo() << "OpenMP worker affinity (this process):"
            << parallel::Pinning::maskToString(seissol::parallel::Pinning::getWorkerUnionMask());
  logInfo() << "OpenMP worker affinity (this node)   :"
            << parallel::Pinning::maskToString(seissol::parallel::Pinning::getNodeMask());

  seissol::printCommThreadInfo(seissol::Mpi::mpi, env_);
  if (seissol::useCommThread(seissol::Mpi::mpi, env_)) {
    auto freeCpus = pinning_.getFreeCPUsMask();
    logInfo() << "Communication thread affinity        :"
              << parallel::Pinning::maskToString(freeCpus);
    if (parallel::Pinning::freeCPUsMaskEmpty(freeCpus)) {
      logError()
          << "There are no free CPUs left. Make sure to leave one for the communication thread. If "
             "you want to run SeisSol without a communication thread (and instead use polling), "
             "then try running with the environment variable \"SEISSOL_COMMTHREAD=0\". ";
    }
  }

  // Check if the ulimit for the stacksize is reasonable.
  // A low limit can lead to segmentation faults.
  rlimit rlim{};
  if (getrlimit(RLIMIT_STACK, &rlim) == 0) {
    const auto rlimInKb = rlim.rlim_cur / 1024;
    // Softlimit (rlim_cur) is enforced by the kernel.
    // This limit is pretty arbitrarily set to 2GiB.
    constexpr auto ReasonableStackLimitInKb = 0x200'000ULL;                    // [kiB] (2 GiB)
    constexpr auto ReasonableStackLimit = ReasonableStackLimitInKb * 0x400ULL; // [B] (2 GiB)
    if (rlim.rlim_cur == RLIM_INFINITY) {
      logInfo() << "The stack size ulimit is unlimited.";
    } else {
      logInfo() << "The stack size ulimit is" << rlimInKb
                << "[kiB] ( =" << UnitByte.formatPrefix(rlim.rlim_cur).c_str() << ").";
    }
    if (rlim.rlim_cur < ReasonableStackLimit) {
      logWarning()
          << "Stack size of" << rlimInKb
          << "[kiB] ( =" << UnitByte.formatPrefix(rlim.rlim_cur).c_str()
          << ") is lower than recommended minimum of" << ReasonableStackLimitInKb
          << "[kiB] ( =" << UnitByte.formatPrefix(ReasonableStackLimit).c_str() << ")."
          << "You can increase the stack size by running the command: ulimit -Ss unlimited.";
    }
  } else {
    logError() << "Stack size cannot be determined because getrlimit syscall failed!";
  }

  // Call post MPI initialization hooks
  seissol::Modules::callHook<ModuleHook::PostMPIInit>();

  // Initialize the ASYNC I/O library
  if (!asyncIO_.init()) {
    return false;
  }

  memoryManager_->initialize();

  return true;
}

void SeisSol::finalize() {
  // Cleanup ASYNC I/O library
  asyncIO_.finalize();

  Modules::callHook<ModuleHook::Shutdown>();

  timeManager_.freeDynamicResources();

  seissol::Mpi::finalize();

  logInfo() << "SeisSol done. Goodbye.";
}

void SeisSol::setBackupTimeStamp(const std::string& stamp) {
  backupTimeStamp_ = stamp;
  seissol::Mpi::mpi.broadcastContainer(backupTimeStamp_, 0);
}

void SeisSol::loadCheckpoint(const std::string& file) {
  checkpointLoadFile_ = std::make_optional<std::string>(file);
}

void SeisSol::setExecutionPlaceCutoff(std::size_t size) { executionPlaceCutoff_ = size; }

} // namespace seissol
