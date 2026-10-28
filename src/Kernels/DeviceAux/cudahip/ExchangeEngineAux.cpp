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

#include <cstring>
#include <hsa/hsa.h>
#include <hsa/hsa_ext_amd.h>
#include <hsa/hsa_ven_amd_loader.h>
#include <string>
#include <utils/logger.h>
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
 * Copies the regions with all threads of the block (of `Threads` threads), in words of 16 bytes
 * where both sides are aligned to them; past the caches that may hold what the source had before
 * (written by other kernels, or by the host).
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
    for (std::size_t i = threadIdx.x; i < words; i += Threads) {
      const auto* from = reinterpret_cast<const uint4*>(source) + i;
#ifdef __HIP__
      // (host memory is not cached; the caches of the block were invalidated before)
      const uint4 value = *from;
#else
      const uint4 value = fromHost ? __ldcv(from) : __ldcg(from);
#endif
      reinterpret_cast<uint4*>(target)[i] = value;
    }
    for (std::size_t i = words * sizeof(uint4) + threadIdx.x; i < region.bytes; i += Threads) {
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

} // namespace

/// the arguments of the engine, in one piece (for its dispatch with HSA)
struct EngineArguments {
  Transport* transports{nullptr};
  const std::uint32_t* control{nullptr};
  std::uint32_t count{0};
};

/**
 * The engine (see run()): with C linkage, so that it can be found by its name among the kernels
 * that HIP has loaded; and with blocks of `Threads` threads, without reading their size (its
 * dispatch with HSA leaves out the hidden arguments of HIP).
 */
extern "C" __global__ void seissolExchangeEngine(EngineArguments arguments) {
  auto* transports = arguments.transports;
  const auto* control = arguments.control;
  const auto count = arguments.count;
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

namespace {
#ifdef __HIP__
/**
 * An HSA queue of its own for the engine: CLR puts all of its streams on a few queues, of whatever
 * priority, so that a kernel that keeps running on a stream holds up the streams that share its
 * queue. The engine gets dispatched there by hand, with the kernel object that HIP has loaded for
 * it; and the synchronization of the device does not wait for it.
 */
class EngineQueue {
  public:
  static EngineQueue& instance() {
    static EngineQueue queue;
    return queue;
  }

  void dispatch(const EngineArguments& arguments, unsigned blocks) {
    if (queue_ == nullptr) {
      setUp();
    }
    wait();
    hsa_signal_store_relaxed(done_, 1);
    std::memcpy(kernargs_, &arguments, sizeof(arguments));
    const auto index = hsa_queue_add_write_index_relaxed(queue_, 1);
    auto* packet =
        static_cast<hsa_kernel_dispatch_packet_t*>(queue_->base_address) + (index % queue_->size);
    // all but the header, which makes the packet valid
    std::memset(reinterpret_cast<std::uint8_t*>(packet) + sizeof(std::uint32_t),
                0,
                sizeof(*packet) - sizeof(std::uint32_t));
    packet->workgroup_size_x = Threads;
    packet->workgroup_size_y = 1;
    packet->workgroup_size_z = 1;
    packet->grid_size_x = Threads * blocks;
    packet->grid_size_y = 1;
    packet->grid_size_z = 1;
    packet->private_segment_size = privateSize_;
    packet->group_segment_size = groupSize_;
    packet->kernel_object = kernelObject_;
    packet->kernarg_address = kernargs_;
    packet->completion_signal = done_;
    const std::uint32_t header =
        (HSA_PACKET_TYPE_KERNEL_DISPATCH << HSA_PACKET_HEADER_TYPE) |
        (HSA_FENCE_SCOPE_SYSTEM << HSA_PACKET_HEADER_SCACQUIRE_FENCE_SCOPE) |
        (HSA_FENCE_SCOPE_SYSTEM << HSA_PACKET_HEADER_SCRELEASE_FENCE_SCOPE);
    const std::uint32_t setup = 1U << HSA_KERNEL_DISPATCH_PACKET_SETUP_DIMENSIONS;
    __atomic_store_n(
        reinterpret_cast<std::uint32_t*>(packet), header | (setup << 16U), __ATOMIC_RELEASE);
    hsa_signal_store_screlease(queue_->doorbell_signal, static_cast<hsa_signal_value_t>(index));
  }

  void wait() {
    if (queue_ != nullptr) {
      while (hsa_signal_wait_scacquire(
                 done_, HSA_SIGNAL_CONDITION_LT, 1, UINT64_MAX, HSA_WAIT_STATE_BLOCKED) >= 1) {
      }
    }
  }

  private:
  static void check(hsa_status_t status, const char* call) {
    if (status != HSA_STATUS_SUCCESS) {
      const char* text = nullptr;
      hsa_status_string(status, &text);
      logError() << "The queue of the exchange engine: HSA failed in" << call << "with"
                 << (text != nullptr ? text : "an unknown error");
    }
  }

  static hsa_status_t visitAgent(hsa_agent_t agent, void* data) {
    auto& queue = *static_cast<EngineQueue*>(data);
    hsa_device_type_t type{};
    hsa_agent_get_info(agent, HSA_AGENT_INFO_DEVICE, &type);
    if (type == HSA_DEVICE_TYPE_CPU && queue.cpu_.handle == 0) {
      queue.cpu_ = agent;
    }
    if (type == HSA_DEVICE_TYPE_GPU) {
      // the device of HIP, by its place on the PCI bus
      std::uint32_t location = 0;
      std::uint32_t domain = 0;
      hsa_agent_get_info(agent, static_cast<hsa_agent_info_t>(HSA_AMD_AGENT_INFO_BDFID), &location);
      hsa_agent_get_info(agent, static_cast<hsa_agent_info_t>(HSA_AMD_AGENT_INFO_DOMAIN), &domain);
      if (((location >> 8U) & 0xffU) == queue.bus_ && ((location >> 3U) & 0x1fU) == queue.slot_ &&
          domain == queue.domain_) {
        queue.gpu_ = agent;
      }
    }
    return HSA_STATUS_SUCCESS;
  }

  static hsa_status_t visitPool(hsa_amd_memory_pool_t pool, void* data) {
    auto& queue = *static_cast<EngineQueue*>(data);
    std::uint32_t flags = 0;
    hsa_amd_memory_pool_get_info(pool, HSA_AMD_MEMORY_POOL_INFO_GLOBAL_FLAGS, &flags);
    if ((flags & HSA_AMD_MEMORY_POOL_GLOBAL_FLAG_KERNARG_INIT) != 0 && !queue.havePool_) {
      queue.kernargPool_ = pool;
      queue.havePool_ = true;
    }
    return HSA_STATUS_SUCCESS;
  }

  static hsa_status_t visitSymbol(hsa_executable_t /*executable*/,
                                  hsa_agent_t /*agent*/,
                                  hsa_executable_symbol_t symbol,
                                  void* data) {
    auto& queue = *static_cast<EngineQueue*>(data);
    hsa_symbol_kind_t kind{};
    hsa_executable_symbol_get_info(symbol, HSA_EXECUTABLE_SYMBOL_INFO_TYPE, &kind);
    if (kind != HSA_SYMBOL_KIND_KERNEL) {
      return HSA_STATUS_SUCCESS;
    }
    std::uint32_t length = 0;
    hsa_executable_symbol_get_info(symbol, HSA_EXECUTABLE_SYMBOL_INFO_NAME_LENGTH, &length);
    std::string name(length, '\0');
    hsa_executable_symbol_get_info(symbol, HSA_EXECUTABLE_SYMBOL_INFO_NAME, name.data());
    // (the name of the kernel descriptor, with the suffix .kd)
    if (name == "seissolExchangeEngine" || name == "seissolExchangeEngine.kd") {
      hsa_executable_symbol_get_info(
          symbol, HSA_EXECUTABLE_SYMBOL_INFO_KERNEL_OBJECT, &queue.kernelObject_);
      hsa_executable_symbol_get_info(
          symbol, HSA_EXECUTABLE_SYMBOL_INFO_KERNEL_KERNARG_SEGMENT_SIZE, &queue.kernargSize_);
      hsa_executable_symbol_get_info(
          symbol, HSA_EXECUTABLE_SYMBOL_INFO_KERNEL_GROUP_SEGMENT_SIZE, &queue.groupSize_);
      hsa_executable_symbol_get_info(
          symbol, HSA_EXECUTABLE_SYMBOL_INFO_KERNEL_PRIVATE_SEGMENT_SIZE, &queue.privateSize_);
    }
    return HSA_STATUS_SUCCESS;
  }

  static hsa_status_t visitExecutable(hsa_executable_t executable, void* data) {
    auto& queue = *static_cast<EngineQueue*>(data);
    hsa_executable_iterate_agent_symbols(executable, queue.gpu_, visitSymbol, data);
    return HSA_STATUS_SUCCESS;
  }

  void setUp() {
    // HIP loads the code object of the engine (its loading may be deferred until it is needed)
    hipFuncAttributes attributes{};
    if (hipFuncGetAttributes(&attributes, reinterpret_cast<const void*>(seissolExchangeEngine)) !=
        hipSuccess) {
      logError() << "The exchange engine could not be loaded.";
    }
    int device = 0;
    hipGetDevice(&device);
    int bus = 0;
    int slot = 0;
    int domain = 0;
    hipDeviceGetAttribute(&bus, hipDeviceAttributePciBusId, device);
    hipDeviceGetAttribute(&slot, hipDeviceAttributePciDeviceId, device);
    hipDeviceGetAttribute(&domain, hipDeviceAttributePciDomainID, device);
    bus_ = static_cast<std::uint32_t>(bus);
    slot_ = static_cast<std::uint32_t>(slot);
    domain_ = static_cast<std::uint32_t>(domain);
    // (HIP has initialized HSA already; this only counts up)
    check(hsa_init(), "hsa_init");
    check(hsa_iterate_agents(visitAgent, this), "hsa_iterate_agents");
    if (gpu_.handle == 0 || cpu_.handle == 0) {
      logError() << "The exchange engine found no HSA agent for its device.";
    }
    check(hsa_amd_agent_iterate_memory_pools(cpu_, visitPool, this),
          "hsa_amd_agent_iterate_memory_pools");
    hsa_ven_amd_loader_1_03_pfn_t loader{};
    check(
        hsa_system_get_major_extension_table(HSA_EXTENSION_AMD_LOADER, 1, sizeof(loader), &loader),
        "hsa_system_get_major_extension_table");
    check(loader.hsa_ven_amd_loader_iterate_executables(visitExecutable, this),
          "hsa_ven_amd_loader_iterate_executables");
    if (kernelObject_ == 0 || !havePool_ || kernargSize_ < sizeof(EngineArguments)) {
      logError() << "The exchange engine found no kernel object of its own (or no memory for its"
                 << "arguments).";
    }
    check(hsa_amd_memory_pool_allocate(kernargPool_, kernargSize_, 0, &kernargs_),
          "hsa_amd_memory_pool_allocate");
    check(hsa_amd_agents_allow_access(1, &gpu_, nullptr, kernargs_), "hsa_amd_agents_allow_access");
    std::memset(kernargs_, 0, kernargSize_);
    check(hsa_signal_create(0, 0, nullptr, &done_), "hsa_signal_create");
    check(hsa_queue_create(gpu_,
                           QueueSize,
                           HSA_QUEUE_TYPE_SINGLE,
                           nullptr,
                           nullptr,
                           UINT32_MAX,
                           UINT32_MAX,
                           &queue_),
          "hsa_queue_create");
  }

  // (one dispatch at a time)
  static constexpr std::uint32_t QueueSize = 64;

  std::uint32_t bus_{0};
  std::uint32_t slot_{0};
  std::uint32_t domain_{0};
  hsa_agent_t gpu_{};
  hsa_agent_t cpu_{};
  hsa_amd_memory_pool_t kernargPool_{};
  bool havePool_{false};
  std::uint64_t kernelObject_{0};
  std::uint32_t kernargSize_{0};
  std::uint32_t groupSize_{0};
  std::uint32_t privateSize_{0};
  hsa_queue_t* queue_{nullptr};
  void* kernargs_{nullptr};
  hsa_signal_t done_{};
};
#endif
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
  const EngineArguments arguments{transports, control, count};
#ifdef __HIP__
  EngineQueue::instance().dispatch(arguments, blocks);
#else
  seissolExchangeEngine<<<blocks, Threads, 0, reinterpret_cast<StreamT>(stream)>>>(arguments);
#endif
}

void waitForStop(void* stream) {
#ifdef __HIP__
  EngineQueue::instance().wait();
#else
  cudaStreamSynchronize(reinterpret_cast<StreamT>(stream));
#endif
}

} // namespace seissol::kernels::device::aux::exchangeengine
