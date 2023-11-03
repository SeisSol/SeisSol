#ifndef SEISSOL_PARALLEL_HELPER_HPP_
#define SEISSOL_PARALLEL_HELPER_HPP_

#include "utils/env.h"

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
} // namespace seissol

#endif // SEISSOL_PARALLEL_HELPER_HPP_
