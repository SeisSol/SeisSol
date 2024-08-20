#include "Parallel/MPI.h"

#ifdef USE_CCL
#ifdef CCL_NCCL
#include <nccl.h>
#endif
#ifdef CCL_RCCL
#include <rccl/rccl.h>
#endif
#include <device.h>
#endif

namespace seissol::solver::clustering::communication {

std::vector<void*> createComms(std::size_t count) {
#ifdef USE_CCL
  // cf. partially https://docs.nvidia.com/deeplearning/nccl/user-guide/docs/examples.html
  std::vector<ncclUniqueId> cclIds(count);
  std::vector<void*> comms(count);
  if (seissol::MPI::mpi.rank() == 0) {
    for (std::size_t i = 0; i < count; ++i) {
      ncclGetUniqueId(&cclIds[i]);
    }
  }
  MPI_Bcast(
      cclIds.data(), sizeof(ncclUniqueId) * cclIds.size(), MPI_BYTE, 0, seissol::MPI::mpi.comm());
  MPI_Barrier(seissol::MPI::mpi.comm());
  for (std::size_t i = 0; i < count; ++i) {
    ncclComm_t preComm;
    ncclCommInitRank(&preComm, seissol::MPI::mpi.size(), cclIds[i], seissol::MPI::mpi.rank());
    comms[i] = static_cast<void*>(preComm);
  }
  return comms;
#else
  return {};
#endif
}

} // namespace seissol::solver::clustering::communication
