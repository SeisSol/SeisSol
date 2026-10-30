// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "RocshmemKernels.h"

#include <cstddef>
#include <cstdint>
#include <hip/hip_runtime.h>
#include <mpi.h>
#include <rocshmem/rocshmem.hpp>

namespace seissol::solver::rocshmem {

namespace {

hipStream_t native(void* stream) { return static_cast<hipStream_t>(stream); }

// a whole work-group puts the data of a region
constexpr unsigned PutThreads = 256;

// the released versions wait on no uint64_t; unsigned long long has the same width
static_assert(sizeof(unsigned long long) == sizeof(std::uint64_t));

__global__ void advanceKernel(std::uint64_t* count) { *count += 1; }

__global__ void signalCountKernel(std::uint64_t* signal, const std::uint64_t* count, int pe) {
  ::rocshmem::rocshmem_uint64_atomic_set(signal, *count, pe);
}

__global__ void
    waitCountKernel(std::uint64_t* signal, const std::uint64_t* count, std::uint64_t factor) {
  ::rocshmem::rocshmem_ulonglong_wait_until(
      reinterpret_cast<unsigned long long*>(signal), ::rocshmem::ROCSHMEM_CMP_GE, *count * factor);
}

__global__ void putSignalKernel(
    void* destination, const void* source, std::size_t bytes, std::uint64_t* signal, int pe) {
  ::rocshmem::rocshmem_putmem_signal_nbi_wg(
      destination, source, bytes, signal, 1, ::rocshmem::ROCSHMEM_SIGNAL_ADD, pe);
}

__global__ void quietKernel() { ::rocshmem::rocshmem_quiet(); }

// copies with the threads of the workgroup, by 16 bytes where it can
__device__ void copyInGroup(void* destination, const void* source, std::uint64_t bytes) {
  const auto start = static_cast<std::uint64_t>(threadIdx.x);
  const auto stride = static_cast<std::uint64_t>(blockDim.x);
  if (((reinterpret_cast<std::uintptr_t>(destination) | reinterpret_cast<std::uintptr_t>(source) |
        bytes) &
       15U) == 0) {
    auto* to = static_cast<uint4*>(destination);
    const auto* from = static_cast<const uint4*>(source);
    for (std::uint64_t i = start; i < bytes / 16; i += stride) {
      to[i] = from[i];
    }
  } else {
    auto* to = static_cast<char*>(destination);
    const auto* from = static_cast<const char*>(source);
    for (std::uint64_t i = start; i < bytes; i += stride) {
      to[i] = from[i];
    }
  }
}

// a group: each workgroup one task (see ShmemGroup.h)
__global__ void groupKernel(const shmem::Task* tasks,
                            const shmem::Region* regions,
                            std::uint64_t* count,
                            std::uint32_t* finished) {
  __shared__ std::uint64_t current;
  if (threadIdx.x == 0) {
    current = __hip_atomic_load(count, __ATOMIC_RELAXED, __HIP_MEMORY_SCOPE_AGENT) + 1;
  }
  __syncthreads();
  const auto& task = tasks[blockIdx.x];
  if (task.receive != 0) {
    if (threadIdx.x == 0) {
      // the staging window is free again: the group before has copied its data out
      ::rocshmem::rocshmem_uint64_atomic_set(task.clearToSend, current, task.peer);
      ::rocshmem::rocshmem_ulonglong_wait_until(reinterpret_cast<unsigned long long*>(task.arrived),
                                                ::rocshmem::ROCSHMEM_CMP_GE,
                                                current * task.regionCount);
    }
    __syncthreads();
    // (the data came from other processes)
    __builtin_amdgcn_fence(__ATOMIC_ACQUIRE, "");
    for (std::uint32_t i = 0; i < task.regionCount; ++i) {
      const auto& region = regions[task.firstRegion + i];
      copyInGroup(region.destination, region.source, region.bytes);
    }
  } else {
    if (threadIdx.x == 0) {
      ::rocshmem::rocshmem_ulonglong_wait_until(
          reinterpret_cast<unsigned long long*>(task.clearToSend),
          ::rocshmem::ROCSHMEM_CMP_GE,
          current);
    }
    __syncthreads();
    for (std::uint32_t i = 0; i < task.regionCount; ++i) {
      const auto& region = regions[task.firstRegion + i];
      ::rocshmem::rocshmem_putmem_signal_nbi_wg(region.destination,
                                                region.source,
                                                region.bytes,
                                                task.arrived,
                                                1,
                                                ::rocshmem::ROCSHMEM_SIGNAL_ADD,
                                                task.peer);
    }
    __syncthreads();
    if (threadIdx.x == 0) {
      // the copy layer may overwrite its data once the puts have left
      ::rocshmem::rocshmem_quiet();
    }
  }
  // the last workgroup advances the count; all of them have read it by then
  __syncthreads();
  if (threadIdx.x == 0) {
    __threadfence();
    if (atomicAdd(finished, 1U) + 1 == gridDim.x) {
      *finished = 0;
      __hip_atomic_store(count, current, __ATOMIC_RELEASE, __HIP_MEMORY_SCOPE_AGENT);
    }
  }
}

} // namespace

int initialize(MPI_Comm comm) {
  ::rocshmem::rocshmem_init_attr_t attributes{};
#ifdef HAVE_EXTERNAL_MPI
  // (read during the initialization)
  attributes.mpi_comm = static_cast<void*>(&comm);
  const auto flags = ::rocshmem::ROCSHMEM_INIT_WITH_MPI_COMM;
#else
  // Built without an MPI at hand (as in the nightlies of TheRock in October 2026, version 3.8),
  // ROCSHMEM looks up the predefined objects of Open MPI with dlsym in the MPI library. That misses
  // the copies an executable gets by copy relocations where it uses these objects itself, as
  // SeisSol does, which are the ones MPI sets up: the check of the IPC backend crashed in
  // MPI_Comm_split_type with MPI_INFO_NULL. With a unique ID that the processes share here,
  // ROCSHMEM bootstraps over TCP instead and leaves these objects alone.
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  ::rocshmem::rocshmem_uniqueid_t uid{};
  if (rank == 0 && ::rocshmem::rocshmem_get_uniqueid(&uid) != ::rocshmem::ROCSHMEM_SUCCESS) {
    return 1;
  }
  MPI_Bcast(uid.data(), static_cast<int>(uid.size()), MPI_BYTE, 0, comm);
  ::rocshmem::rocshmem_set_attr_uniqueid_args(rank, size, &uid, &attributes);
  const auto flags = ::rocshmem::ROCSHMEM_INIT_WITH_UNIQUEID;
#endif
  return ::rocshmem::rocshmem_init_attr(flags, &attributes) == ::rocshmem::ROCSHMEM_SUCCESS ? 0 : 1;
}

void finalize() { ::rocshmem::rocshmem_finalize(); }

int pe() { return ::rocshmem::rocshmem_my_pe(); }

void* allocate(std::size_t bytes) { return ::rocshmem::rocshmem_malloc(bytes); }

void deallocate(void* pointer) { ::rocshmem::rocshmem_free(pointer); }

void barrierAll() { ::rocshmem::rocshmem_barrier_all(); }

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
           std::uint64_t* count,
           std::uint32_t* finished,
           void* stream) {
  groupKernel<<<static_cast<unsigned>(taskCount), PutThreads, 0, native(stream)>>>(
      tasks, regions, count, finished);
}

const char* lastError() {
  const auto error = hipGetLastError();
  return error == hipSuccess ? nullptr : hipGetErrorString(error);
}

} // namespace seissol::solver::rocshmem
