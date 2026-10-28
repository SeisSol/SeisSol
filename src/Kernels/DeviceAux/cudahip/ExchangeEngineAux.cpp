// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/DeviceAux/ExchangeEngineAux.h"

#include <cstddef>
#include <cstdint>

#ifdef __HIP__
#include "hip/hip_runtime.h"
using StreamT = hipStream_t;
#endif
#ifdef __CUDACC__
using StreamT = cudaStream_t;
#endif

namespace seissol::kernels::device::aux::exchangeengine {

namespace {
constexpr unsigned Threads = 256;
// the longest pause between two looks for work: about 5 us
constexpr unsigned MaxPause = 128;

__device__ std::uint32_t loadDevice(const std::uint32_t* value) {
#ifdef __HIP__
  return __hip_atomic_load(value, __ATOMIC_ACQUIRE, __HIP_MEMORY_SCOPE_AGENT);
#else
  const auto loaded = *static_cast<const volatile std::uint32_t*>(value);
  __threadfence();
  return loaded;
#endif
}

__device__ std::uint32_t loadHost(const std::uint32_t* value) {
#ifdef __HIP__
  return __hip_atomic_load(value, __ATOMIC_ACQUIRE, __HIP_MEMORY_SCOPE_SYSTEM);
#else
  const auto loaded = *static_cast<const volatile std::uint32_t*>(value);
  __threadfence_system();
  return loaded;
#endif
}

__device__ void storeHost(std::uint32_t* target, std::uint32_t value) {
  // what the block has written before reaches the host first
  __threadfence_system();
#ifdef __HIP__
  __hip_atomic_store(target, value, __ATOMIC_RELEASE, __HIP_MEMORY_SCOPE_SYSTEM);
#else
  *static_cast<volatile std::uint32_t*>(target) = value;
#endif
}

__device__ void pauseFor(unsigned pause) {
#ifdef __HIP__
  // about 64 clocks each
  for (unsigned i = 0; i < pause; ++i) {
    __builtin_amdgcn_s_sleep(1);
  }
#else
  __nanosleep(32 * pause);
#endif
}

/**
 * Copies the regions with all threads of the block, in words of 16 bytes where both sides are
 * aligned to them; past the caches that may hold what the source had before (written by other
 * kernels, or by the host).
 */
__device__ void copyRegions(const Region* regions, std::uint32_t count, bool fromHost) {
  for (std::uint32_t current = 0; current < count; ++current) {
    const auto region = regions[current];
    auto* target = static_cast<unsigned char*>(region.target);
    const auto* source = static_cast<const unsigned char*>(region.source);
    const bool aligned =
        (reinterpret_cast<std::uintptr_t>(target) | reinterpret_cast<std::uintptr_t>(source)) %
            sizeof(uint4) ==
        0;
    const std::size_t words = aligned ? region.bytes / sizeof(uint4) : 0;
    for (std::size_t i = threadIdx.x; i < words; i += blockDim.x) {
      const auto* from = reinterpret_cast<const uint4*>(source) + i;
#ifdef __HIP__
      // (host memory is not cached; the caches of the block were invalidated before)
      const uint4 value = *from;
#else
      const uint4 value = fromHost ? __ldcv(from) : __ldcg(from);
#endif
      reinterpret_cast<uint4*>(target)[i] = value;
    }
    for (std::size_t i = words * sizeof(uint4) + threadIdx.x; i < region.bytes; i += blockDim.x) {
      target[i] = static_cast<const volatile unsigned char*>(source)[i];
    }
  }
}

__global__ void requestKernel(Transport* sender,
                              std::uint32_t send,
                              Transport* receiver,
                              std::uint32_t receive) {
  // (the work of the stream before has completed, and is there for the engine)
  __threadfence();
  if (sender != nullptr) {
    atomicMax(&sender->sendRequested, send);
  }
  if (receiver != nullptr) {
    atomicMax(&receiver->receiveRequested, receive);
  }
}

__global__ void
    runKernel(Transport* transports, std::uint32_t count, const std::uint32_t* control) {
  // the send (2 t) or receive (2 t + 1) of transport t that the block takes up, -1 for none; and
  // its number
  __shared__ int job;
  __shared__ std::uint32_t number;
  __shared__ bool stop;
  unsigned pause = 1;
  // where the block starts to look, so that the blocks do not all look at the same transports
  std::uint32_t start = blockIdx.x % count;
  while (true) {
    if (threadIdx.x == 0) {
      job = -1;
      stop = false;
      for (std::uint32_t i = 0; i < count; ++i) {
        const auto index = (start + i) % count;
        auto& transport = transports[index];
        // a send once requested, and once MPI is done with the one before
        const auto send = loadDevice(&transport.sendClaimed) + 1;
        if (loadDevice(&transport.sendRequested) >= send && loadHost(transport.sent) + 1 >= send &&
            atomicCAS(&transport.sendClaimed, send - 1, send) == send - 1) {
          job = static_cast<int>(2 * index);
          number = send;
          break;
        }
        // a receive once requested, and once MPI has received it
        const auto receive = loadDevice(&transport.receiveClaimed) + 1;
        if (loadDevice(&transport.receiveRequested) >= receive &&
            loadHost(transport.received) >= receive &&
            atomicCAS(&transport.receiveClaimed, receive - 1, receive) == receive - 1) {
          job = static_cast<int>(2 * index + 1);
          number = receive;
          break;
        }
      }
      start = (start + 1) % count;
      if (job < 0 && loadHost(control) != 0) {
        bool done = true;
        for (std::uint32_t index = 0; index < count && done; ++index) {
          done = loadDevice(&transports[index].sendClaimed) >= loadHost(control + 1 + index) &&
                 loadDevice(&transports[index].receiveClaimed) >=
                     loadHost(control + 1 + count + index);
        }
        stop = done;
      }
    }
    __syncthreads();
    if (stop) {
      return;
    }
    if (job < 0) {
      if (threadIdx.x == 0) {
        pauseFor(pause);
      }
      pause = pause < MaxPause ? 2 * pause : MaxPause;
      __syncthreads();
      continue;
    }
    pause = 1;
#ifdef __HIP__
    // what the source had before is not in the caches of the block any more
    __builtin_amdgcn_fence(__ATOMIC_ACQUIRE, "agent");
#endif
    auto& transport = transports[job / 2];
    const bool send = job % 2 == 0;
    if (send) {
      copyRegions(transport.out, transport.outCount, false);
    } else {
      copyRegions(transport.in, transport.inCount, true);
    }
    __syncthreads();
    if (threadIdx.x == 0) {
      storeHost(send ? transport.copiedOut : transport.copiedIn, number);
    }
    __syncthreads();
  }
}
} // namespace

void request(Transport* sender,
             std::uint32_t send,
             Transport* receiver,
             std::uint32_t receive,
             void* stream) {
  if (sender == nullptr && receiver == nullptr) {
    return;
  }
  requestKernel<<<1, 1, 0, reinterpret_cast<StreamT>(stream)>>>(sender, send, receiver, receive);
}

void run(Transport* transports,
         std::uint32_t count,
         const std::uint32_t* control,
         unsigned blocks,
         void* stream) {
  if (count == 0) {
    return;
  }
  runKernel<<<blocks, Threads, 0, reinterpret_cast<StreamT>(stream)>>>(transports, count, control);
}

} // namespace seissol::kernels::device::aux::exchangeengine
