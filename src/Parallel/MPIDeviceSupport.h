// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_MPIDEVICESUPPORT_H_
#define SEISSOL_SRC_PARALLEL_MPIDEVICESUPPORT_H_

#include "Parallel/MPI.h"

#include <optional>
#include <string>

namespace seissol::parallel {

/**
 * What the MPI library tells about whether it can read and write the device memory of this build.
 */
struct MpiDeviceSupport {
  /// whether it can; none if it does not tell, e.g. without a query for it, or in a SYCL build,
  /// which does not know the device memory of which backend it works on
  std::optional<bool> supported;
  /// what has told it, e.g. "MPIX_Query_cuda_support"; empty if nothing has
  std::string source;
};

/**
 * Asks the MPI library of this process (after MPI_Init): Open MPI with MPIX_Query_cuda_support or
 * MPIX_Query_rocm_support (or by MPIX_CUDA_AWARE_SUPPORT, MPIX_ROCM_AWARE_SUPPORT, if it has been
 * built without), MPICH and the libraries that derive from it (e.g. Intel MPI) with
 * MPIX_GPU_query_support.
 */
MpiDeviceSupport queryMpiDeviceSupport();

/**
 * The transfer mode of the halo exchange: the one asked for (SEISSOL_TRANSFER_MODE), if any;
 * otherwise direct, unless the MPI library tells that it cannot read and write device memory, and
 * the buffers of the exchange are device memory (i.e. not unified memory, SEISSOL_USM_MPI): then
 * host, which stages the exchange in host memory.
 */
Mpi::DataTransferMode chooseDataTransferMode(std::optional<Mpi::DataTransferMode> requested,
                                             std::optional<bool> deviceSupport,
                                             bool unifiedBuffers);

/**
 * Chooses the transfer mode as chooseDataTransferMode() does, with the support of the MPI library
 * as the ranks of `comm` have found it (the least of them, so that all choose the same), and
 * reports what it has found; collective.
 */
Mpi::DataTransferMode decideDataTransferMode(std::optional<Mpi::DataTransferMode> requested,
                                             bool unifiedBuffers,
                                             bool sharedDeviceMemory,
                                             MPI_Comm comm);

} // namespace seissol::parallel

#endif // SEISSOL_SRC_PARALLEL_MPIDEVICESUPPORT_H_
