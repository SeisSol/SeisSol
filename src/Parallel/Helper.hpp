#ifndef SEISSOL_PARALLEL_HELPER_HPP_
#define SEISSOL_PARALLEL_HELPER_HPP_

#include "utils/env.h"

#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol {
template <typename T>
void printCommThreadInfo(const T& mpiBasic) {
  bool useThread = utils::Env::get<bool>("SEISSOL_COMMTHREAD", true);
  if (mpiBasic.isSingleProcess()) {
    logInfo(mpiBasic.rank()) << "Using polling for advancing MPI communication, due to having only "
                                "one MPI rank running.";
  } else if (useThread) {
    logInfo(mpiBasic.rank()) << "Using a communication thread for advancing MPI communication.";
  } else {
    logInfo(mpiBasic.rank()) << "Using polling for advancing MPI communication.";
  }
}

template <typename T>
bool useCommThread(const T& mpiBasic) {
  bool useThread = utils::Env::get<bool>("SEISSOL_COMMTHREAD", true);
  return useThread && !mpiBasic.isSingleProcess();
}

inline bool usePersistentMpi() { return utils::Env::get<bool>("SEISSOL_MPI_PERSISTENT", false); }

template <typename T>
void printPersistentMpiInfo(const T& mpiBasic) {
  if (usePersistentMpi()) {
    logInfo(mpiBasic.rank()) << "Using persistent MPI routines.";
  } else {
    logInfo(mpiBasic.rank()) << "Using asynchronous MPI routines.";
  }
}

#ifdef ACL_DEVICE
inline bool useUSM() { return utils::Env::get<bool>("SEISSOL_USM", device::DeviceInstance::getInstance().api->isUnifiedMemoryDefault()); }

template <typename T>
void printUSMInfo(const T& mpiBasic) {
  if (useUSM()) {
    logInfo(mpiBasic.rank()) << "Using unified buffers for CPU-GPU data.";
  } else {
    logInfo(mpiBasic.rank()) << "Using separate buffers for CPU-GPU data.";
  }
}

inline bool useMPIUSM() { return utils::Env::get<bool>("SEISSOL_USM_MPI", device::DeviceInstance::getInstance().api->isUnifiedMemoryDefault()); }

template <typename T>
void printMPIUSMInfo(const T& mpiBasic) {
  if (useMPIUSM()) {
    logInfo(mpiBasic.rank()) << "Using unified buffers for CPU-GPU MPI data.";
  } else {
    logInfo(mpiBasic.rank()) << "Using separate buffers for CPU-GPU MPI data.";
  }
}

inline int deviceHostSwitch() { return utils::Env::get<int>("SEISSOL_DEVICE_HOST_SWITCH", 0); }

template <typename T>
void printDeviceHostSwitch(const T& mpiBasic) {
  logInfo(mpiBasic.rank()) << "Running clusters with" << deviceHostSwitch()
                           << "or more cells on the GPU (and on the CPU otherwise)";
}
#endif

} // namespace seissol

#endif // SEISSOL_PARALLEL_HELPER_HPP_
