// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "NvshmemKernels.h"

#include <cstddef>
#include <cstdint>
#include <cuda.h>
#include <cuda_runtime.h>
#include <map>
#include <mpi.h>
#include <nvshmem.h>
#include <nvshmemx.h>
#include <vector>

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

const char* registerSources(const void* const* sources,
                            const std::size_t* bytes,
                            std::size_t count,
                            std::vector<void*>& registered) {
  // the allocations, by their base addresses; the runtime API cannot tell where one starts
  std::map<CUdeviceptr, std::size_t> allocations;
  for (std::size_t i = 0; i < count; ++i) {
    auto address = reinterpret_cast<CUdeviceptr>(sources[i]);
    const auto end = address + bytes[i];
    while (address < end) {
      CUdeviceptr base = 0;
      std::size_t size = 0;
      if (cuMemGetAddressRange(&base, &size, address) != CUDA_SUCCESS || size == 0) {
        return "the memory of the copy layers is no device allocation that CUDA knows";
      }
      allocations[base] = size;
      address = base + size;
    }
  }
  for (const auto& [base, size] : allocations) {
    auto* pointer = reinterpret_cast<void*>(base);
    if (nvshmemx_buffer_register(pointer, size) != 0) {
      return "nvshmemx_buffer_register failed";
    }
    registered.push_back(pointer);
  }
  return nullptr;
}

void unregisterSources(const std::vector<void*>& registered) {
  for (auto* pointer : registered) {
    nvshmemx_buffer_unregister(pointer);
  }
}

} // namespace seissol::solver::nvshmem
