#include "Parallel/MPI.h"

#ifdef USE_CCL
#ifdef SEISSOL_KERNELS_CUDA
#include <nccl.h>
#endif
#ifdef SEISSOL_KERNELS_HIP
#include <rccl/rccl.h>
#endif
#include <device.h>
#endif

namespace seissol::solver::clustering::communication {

void* createComm() {
#ifdef USE_CCL
  // cf. partially https://docs.nvidia.com/deeplearning/nccl/user-guide/docs/examples.html
  ncclComm_t preComm;
  ncclUniqueId cclId;
  if (seissol::MPI::mpi.rank() == 0) {
    ncclGetUniqueId(&cclId);
  }
  MPI_Bcast(&cclId, sizeof(ncclUniqueId), MPI_BYTE, 0, seissol::MPI::mpi.comm());
  ncclCommInitRank(&preComm, seissol::MPI::mpi.size(), cclId, seissol::MPI::mpi.rank());
  return static_cast<void*>(preComm);
#else
  return nullptr;
#endif
}

} // namespace seissol::solver::clustering::communication
