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

static_assert(sizeof(unsigned long long) == sizeof(std::uint64_t));

// copies with the threads of the block, by 16 bytes where it can; the source is what the network
// card wrote, read past the cache of the SM
__device__ void copyInBlock(void* destination, const void* source, std::uint64_t bytes) {
  const auto start = static_cast<std::uint64_t>(threadIdx.x);
  const auto stride = static_cast<std::uint64_t>(blockDim.x);
  if (((reinterpret_cast<std::uintptr_t>(destination) | reinterpret_cast<std::uintptr_t>(source) |
        bytes) &
       15U) == 0) {
    auto* to = static_cast<uint4*>(destination);
    const auto* from = static_cast<const uint4*>(source);
    for (std::uint64_t i = start; i < bytes / 16; i += stride) {
      to[i] = __ldcg(from + i);
    }
  } else {
    auto* to = static_cast<char*>(destination);
    const auto* from = static_cast<const char*>(source);
    for (std::uint64_t i = start; i < bytes; i += stride) {
      to[i] = __ldcg(from + i);
    }
  }
}

// a group: each block one task (see ShmemGroup.h)
__global__ void groupKernel(const shmem::Task* tasks,
                            const shmem::Region* regions,
                            std::uint64_t windowStride,
                            std::uint64_t* count,
                            std::uint32_t* finished) {
  __shared__ std::uint64_t current;
  if (threadIdx.x == 0) {
    current = atomicAdd(reinterpret_cast<unsigned long long*>(count), 0ULL) + 1;
  }
  __syncthreads();
  const auto& task = tasks[blockIdx.x];
  // the copy of the staging window of this exchange (see ShmemGroup.h)
  const auto window = (current % 2) * windowStride;
  if (task.receive != 0) {
    if (threadIdx.x == 0) {
      nvshmem_signal_wait_until(task.arrived, NVSHMEM_CMP_GE, current * task.regionCount);
      // (the data came from other processes)
      __threadfence_system();
    }
    __syncthreads();
    for (std::uint32_t i = 0; i < task.regionCount; ++i) {
      const auto& region = regions[task.firstRegion + i];
      copyInBlock(
          region.destination, static_cast<const char*>(region.source) + window, region.bytes);
    }
    __syncthreads();
    if (threadIdx.x == 0) {
      // this copy of the window is free again, for the exchange after the next one
      nvshmemx_signal_op(task.clearToSend, current, NVSHMEM_SIGNAL_SET, task.peer);
    }
  } else {
    if (threadIdx.x == 0) {
      // the receiver has copied the exchange before the last one out of this copy of the window
      nvshmem_signal_wait_until(task.clearToSend, NVSHMEM_CMP_GE, current > 2 ? current - 2 : 0);
    }
    __syncthreads();
    for (std::uint32_t i = 0; i < task.regionCount; ++i) {
      const auto& region = regions[task.firstRegion + i];
      nvshmemx_putmem_signal_nbi_block(static_cast<char*>(region.destination) + window,
                                       region.source,
                                       region.bytes,
                                       task.arrived,
                                       1,
                                       NVSHMEM_SIGNAL_ADD,
                                       task.peer);
    }
    __syncthreads();
    if (threadIdx.x == 0) {
      // the copy layer may overwrite its data once the puts have left
      nvshmem_quiet();
    }
  }
  // the last block advances the count; all of them have read it by then
  __syncthreads();
  if (threadIdx.x == 0) {
    __threadfence();
    if (atomicAdd(finished, 1U) + 1 == gridDim.x) {
      *finished = 0;
      __threadfence();
      atomicExch(reinterpret_cast<unsigned long long*>(count), current);
    }
  }
}

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

void group(const shmem::Task* tasks,
           std::size_t taskCount,
           const shmem::Region* regions,
           std::uint64_t windowStride,
           std::uint64_t* count,
           std::uint32_t* finished,
           void* stream) {
  groupKernel<<<static_cast<unsigned>(taskCount), PutThreads, 0, native(stream)>>>(
      tasks, regions, windowStride, count, finished);
}

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
