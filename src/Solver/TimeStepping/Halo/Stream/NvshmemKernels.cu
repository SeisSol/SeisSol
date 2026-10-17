// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "NvshmemKernels.h"

#include <cstddef>
#include <cstdint>
#include <cuda_runtime.h>
#include <mpi.h>
#include <nvshmem.h>
#include <nvshmemx.h>

namespace seissol::solver::nvshmem {

namespace {

cudaStream_t native(void* stream) { return static_cast<cudaStream_t>(stream); }

// a whole block puts the data of a region
constexpr unsigned PutThreads = 512;

__global__ void advanceKernel(std::uint64_t* count) { *count += 1; }

__global__ void signalCountKernel(std::uint64_t* signal, const std::uint64_t* count, int pe) {
  nvshmemx_signal_op(signal, *count, NVSHMEM_SIGNAL_SET, pe);
}

__global__ void
    waitCountKernel(std::uint64_t* signal, const std::uint64_t* count, std::uint64_t factor) {
  nvshmem_signal_wait_until(signal, NVSHMEM_CMP_GE, *count * factor);
}

__global__ void putSignalKernel(
    void* destination, const void* source, std::size_t bytes, std::uint64_t* signal, int pe) {
  nvshmemx_putmem_signal_nbi_block(destination, source, bytes, signal, 1, NVSHMEM_SIGNAL_ADD, pe);
}

__global__ void quietKernel() { nvshmem_quiet(); }

} // namespace

int initialize(MPI_Comm comm) {
  nvshmemx_init_attr_t attributes{};
  nvshmemx_set_attr_mpi_comm_args(&comm, &attributes);
  return nvshmemx_init_attr(NVSHMEMX_INIT_WITH_MPI_COMM, &attributes);
}

void finalize() { nvshmem_finalize(); }

void advance(std::uint64_t* count, void* stream) {
  advanceKernel<<<1, 1, 0, native(stream)>>>(count);
}

void signalCount(std::uint64_t* signal, const std::uint64_t* count, int pe, void* stream) {
  signalCountKernel<<<1, 1, 0, native(stream)>>>(signal, count, pe);
}

void waitCount(std::uint64_t* signal,
               const std::uint64_t* count,
               std::uint64_t factor,
               void* stream) {
  waitCountKernel<<<1, 1, 0, native(stream)>>>(signal, count, factor);
}

void putSignal(void* destination,
               const void* source,
               std::size_t bytes,
               std::uint64_t* signal,
               int pe,
               void* stream) {
  putSignalKernel<<<1, PutThreads, 0, native(stream)>>>(destination, source, bytes, signal, pe);
}

void quiet(void* stream) { quietKernel<<<1, 1, 0, native(stream)>>>(); }

const char* lastError() {
  const auto error = cudaGetLastError();
  return error == cudaSuccess ? nullptr : cudaGetErrorString(error);
}

} // namespace seissol::solver::nvshmem
