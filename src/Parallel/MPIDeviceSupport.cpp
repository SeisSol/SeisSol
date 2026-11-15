// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "MPIDeviceSupport.h"

#include "Parallel/MPI.h"

#include <mpi.h>
#include <optional>
#include <string>
#include <utils/logger.h>

// the extensions of Open MPI, among them its queries for device memory (used in the builds for
// CUDA and HIP only)
#if defined(OPEN_MPI) && __has_include(<mpi-ext.h>)
#include <mpi-ext.h> // IWYU pragma: keep
#endif

namespace seissol::parallel {

MpiDeviceSupport queryMpiDeviceSupport() {
#ifdef SEISSOL_KERNELS_CUDA
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
  return MpiDeviceSupport{.supported = MPIX_Query_cuda_support() != 0,
                          .source = "MPIX_Query_cuda_support"};
#elif defined(MPIX_CUDA_AWARE_SUPPORT)
  // Open MPI, built without CUDA
  return MpiDeviceSupport{.supported = false, .source = "MPIX_CUDA_AWARE_SUPPORT"};
#elif defined(MPIX_GPU_SUPPORT_CUDA)
  int supported = 0;
  MPIX_GPU_query_support(MPIX_GPU_SUPPORT_CUDA, &supported);
  return MpiDeviceSupport{.supported = supported != 0, .source = "MPIX_GPU_query_support"};
#endif
#elif defined(SEISSOL_KERNELS_HIP)
#if defined(MPIX_ROCM_AWARE_SUPPORT) && MPIX_ROCM_AWARE_SUPPORT
  return MpiDeviceSupport{.supported = MPIX_Query_rocm_support() != 0,
                          .source = "MPIX_Query_rocm_support"};
#elif defined(MPIX_ROCM_AWARE_SUPPORT)
  // Open MPI, built without ROCm
  return MpiDeviceSupport{.supported = false, .source = "MPIX_ROCM_AWARE_SUPPORT"};
#elif defined(MPIX_GPU_SUPPORT_HIP)
  int supported = 0;
  MPIX_GPU_query_support(MPIX_GPU_SUPPORT_HIP, &supported);
  return MpiDeviceSupport{.supported = supported != 0, .source = "MPIX_GPU_query_support"};
#endif
#endif
  return {};
}

Mpi::DataTransferMode chooseDataTransferMode(std::optional<Mpi::DataTransferMode> requested,
                                             std::optional<bool> deviceSupport,
                                             bool unifiedBuffers) {
  if (requested.has_value()) {
    return requested.value();
  }
  const bool unsupported = deviceSupport.has_value() && !deviceSupport.value();
  if (unsupported && !unifiedBuffers) {
    return Mpi::DataTransferMode::CopyInCopyOutHost;
  }
  return Mpi::DataTransferMode::Direct;
}

Mpi::DataTransferMode decideDataTransferMode(std::optional<Mpi::DataTransferMode> requested,
                                             bool unifiedBuffers,
                                             bool sharedDeviceMemory,
                                             MPI_Comm comm) {
  const auto support = queryMpiDeviceSupport();
  // the least support any rank has found (none, not told, full), so that all ranks choose the same
  int answer = 1;
  if (support.supported.has_value()) {
    answer = support.supported.value() ? 2 : 0;
  }
  MPI_Allreduce(MPI_IN_PLACE, &answer, 1, MPI_INT, MPI_MIN, comm);
  std::optional<bool> supported;
  if (answer != 1) {
    supported = answer == 2;
  }
  const auto mode = chooseDataTransferMode(requested, supported, unifiedBuffers);

  // (the ranks have the same MPI library, else this rank names what another has found)
  const std::string told =
      support.source.empty() ? std::string("(on another rank)") : "(" + support.source + ")";
  if (!supported.has_value()) {
    if (mode == Mpi::DataTransferMode::Direct && !unifiedBuffers) {
      logInfo() << "The MPI library does not tell whether it can read and write device memory;"
                << "the transfer mode direct assumes that it can (SEISSOL_TRANSFER_MODE=host"
                << "stages the halo exchange in host memory instead).";
    }
  } else if (supported.value()) {
    logInfo() << "The MPI library can read and write device memory" << (told + ".").c_str();
  } else if (unifiedBuffers) {
    logInfo() << "The MPI library cannot read and write device memory" << (told + ";").c_str()
              << "it works on the buffers of the halo exchange nonetheless, as they are"
              << "unified memory (SEISSOL_USM_MPI=1).";
  } else if (mode == Mpi::DataTransferMode::Direct) {
    logWarning() << "The MPI library tells that it cannot read and write device memory"
                 << (told + ",").c_str() << "but gets the device memory of the halo exchange"
                 << "(SEISSOL_TRANSFER_MODE=direct), which may crash it."
                 << "SEISSOL_TRANSFER_MODE=host stages the exchange in host memory instead.";
  } else if (!requested.has_value()) {
    logWarning() << "The MPI library cannot read and write device memory" << (told + ";").c_str()
                 << "hence, the halo exchange goes through host memory (transfer mode host;"
                 << "SEISSOL_TRANSFER_MODE=direct hands it the device memory nonetheless).";
    if (sharedDeviceMemory) {
      logInfo() << "The device shares its memory with the host: with SEISSOL_USM_MPI=1, the MPI"
                << "library works on the buffers of the halo exchange in place instead, without"
                << "copying them to and from host memory.";
    }
  }
  return mode;
}

} // namespace seissol::parallel
