// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#include "MemoryAllocator.h"

#include "Kernels/Common.h"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <stdlib.h>
#include <type_traits>
#include <utils/logger.h>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

#ifdef USE_SHMEM

#ifdef USE_NVSHMEM
#include <nvshmem.h>
#include <nvshmemx.h>
#endif

#ifdef USE_ROCSHMEM
#include <rocshmem.h>
#endif

#ifdef USE_ISHMEM
#include <ishmem.h>
#include <ishmemx.h>
#endif

namespace {

void* mallocShmemInternal(std::size_t size, std::size_t alignment) {
#ifdef USE_NVSHMEM
  return nvshmem_align(alignment, size);
#endif
#ifdef USE_ROCSHMEM
  return rocshmem_align(alignment, size);
#endif
#ifdef USE_ISHMEM
  return ishmem_align(alignment, size);
#endif
}

void* mallocShmem(std::size_t size, std::size_t alignment) {
  std::size_t trueSize = size;

  // shmem_align/shmem_malloc needs to be called with the same size on all ranks
  // (so just take the max of all regions)

  MPI_Allreduce(&size,
                &trueSize,
                1,
                seissol::Mpi::castToMpiType<std::size_t>(),
                MPI_MAX,
                seissol::Mpi::mpi.comm());

  return mallocShmemInternal(trueSize, alignment);
}

void* freeShmemInternal(void* ptr) {
#ifdef USE_NVSHMEM
  nvshmem_free(ptr);
#endif
#ifdef USE_ROCSHMEM
  rocshmem_free(ptr);
#endif
#ifdef USE_ISHMEM
  ishmem_free(ptr);
#endif
}

void freeShmem(void* ptr) { freeShmemInternal(ptr); }

} // namespace

#endif

// TODO: add Device/#48; then refactor away most ifdefs in this file

namespace {
#ifdef USE_MEMKIND
constexpr bool HasHBM = true;
#else
constexpr bool HasHBM = false;
#endif
} // namespace

namespace seissol::memory {

void* allocate(size_t size, size_t alignment, Memkind memkind) {
  // presumably a clang-tidy bug; it'd want to make ptrBuffer const
  // NOLINTNEXTLINE(misc-const-correctness)
  void* ptrBuffer{nullptr};
  bool error = false;

  // handle zero-size allocations
  if (size == 0 && memkind != Memkind::Shmem) {
    ptrBuffer = nullptr;
    return ptrBuffer;
  }

  // switch memkind; reproducing legacy behavior
  if (memkind == Memkind::HighBandwidth && !HasHBM) {
    memkind = Memkind::Standard;
  }
  if ((memkind == Memkind::DeviceUnifiedMemory || memkind == Memkind::PinnedMemory) &&
      !isDeviceOn()) {
    memkind = Memkind::Standard;
  }

  if (memkind == Memkind::Standard) {
    if (alignment % (sizeof(void*)) != 0) {
      ptrBuffer = malloc(size);
      error = (ptrBuffer == nullptr);
    } else {
      error = (posix_memalign(&ptrBuffer, alignment, size) != 0);
    }
  } else if (memkind == Memkind::HighBandwidth) {
#ifdef USE_MEMKIND
    if (alignment % (sizeof(void*)) != 0) {
      ptrBuffer = hbw_malloc(size);
      error = (ptrBuffer == nullptr);
    } else {
      error = (hbw_posix_memalign(&ptrBuffer, alignment, size) != 0);
    }
#endif
  } else if (memkind == Memkind::DeviceGlobalMemory) {
#ifdef ACL_DEVICE
    ptrBuffer = device::DeviceInstance::getInstance().api->allocGlobMem(size);
#endif
  } else if (memkind == Memkind::DeviceUnifiedMemory) {
#ifdef ACL_DEVICE
    ptrBuffer = device::DeviceInstance::getInstance().api->allocUnifiedMem(size);
#endif
  } else if (memkind == Memkind::PinnedMemory) {
#ifdef ACL_DEVICE
    ptrBuffer = device::DeviceInstance::getInstance().api->allocPinnedMem(size);
  } else if (memkind == Memkind::DeviceGlobalCompressed) {
    ptrBuffer = device::DeviceInstance::getInstance().api->allocGlobMem(size, true);
#endif
  } else if (memkind == Memkind::DeviceGlobalCompressed) {
#ifdef ACL_DEVICE
    ptrBuffer = device::DeviceInstance::getInstance().api->allocGlobMem(size, true);
#endif
  } else if (memkind == Memkind::Shmem) {
#ifdef USE_SHMEM
    ptrBuffer = mallocShmem(size);
#endif
  } else {
    logError() << "unknown memkind type used ("
               << static_cast<std::underlying_type_t<Memkind>>(memkind)
               << "). Please, refer to the documentation";
  }

  if (error || ptrBuffer == nullptr) {
    logError() << "The malloc failed (bytes:" << size << ", alignment:" << alignment
               << ", memkind:" << static_cast<std::underlying_type_t<Memkind>>(memkind) << ").";
  }

  return ptrBuffer;
}

void free(void* pointer, Memkind memkind) {

  if (memkind == Memkind::HighBandwidth && !HasHBM) {
    memkind = Memkind::Standard;
  }
  if ((memkind == Memkind::DeviceUnifiedMemory || memkind == Memkind::PinnedMemory) &&
      !isDeviceOn()) {
    memkind = Memkind::Standard;
  }

  if (memkind == Memkind::Standard) {
    ::free(pointer);
  } else if (memkind == Memkind::HighBandwidth) {
#ifdef USE_MEMKIND
    hbw_free(pointer);
#endif
  } else if (memkind == Memkind::DeviceGlobalMemory) {
#ifdef ACL_DEVICE
    device::DeviceInstance::getInstance().api->freeGlobMem(pointer);
#endif
  } else if (memkind == Memkind::DeviceUnifiedMemory) {
#ifdef ACL_DEVICE
    device::DeviceInstance::getInstance().api->freeUnifiedMem(pointer);
#endif
  } else if (memkind == Memkind::PinnedMemory) {
#ifdef ACL_DEVICE
    device::DeviceInstance::getInstance().api->freePinnedMem(pointer);
  } else if (memkind == Memkind::DeviceGlobalCompressed) {
    device::DeviceInstance::getInstance().api->freeGlobMem(pointer);
#endif
  } else if (memkind == Memkind::DeviceGlobalCompressed) {
#ifdef ACL_DEVICE
    device::DeviceInstance::getInstance().api->freeGlobMem(pointer);
#endif
  } else if (memkind == Memkind::Shmem) {
#ifdef USE_SHMEM
    freeShmem(pointer);
#endif
  } else {
    logError() << "unknown memkind type used ("
               << static_cast<std::underlying_type_t<Memkind>>(memkind)
               << "). Please, refer to the documentation";
  }
}

void memcopy(void* dst, const void* src, std::size_t size, Memkind dstMemkind, Memkind srcMemkind) {
  if (dstMemkind == Memkind::DeviceGlobalMemory && srcMemkind != Memkind::DeviceGlobalMemory) {
#ifdef ACL_DEVICE
    device::DeviceInstance::getInstance().api->copyTo(dst, src, size);
#endif
  } else if (dstMemkind != Memkind::DeviceGlobalMemory &&
             srcMemkind == Memkind::DeviceGlobalMemory) {
#ifdef ACL_DEVICE
    device::DeviceInstance::getInstance().api->copyFrom(dst, src, size);
#endif
  } else if (dstMemkind == Memkind::DeviceGlobalMemory &&
             srcMemkind == Memkind::DeviceGlobalMemory) {
#ifdef ACL_DEVICE
    device::DeviceInstance::getInstance().api->copyBetween(dst, src, size);
#endif
  } else {
    std::memcpy(dst, src, size);
  }
}

void memzero(void* dst, std::size_t size, enum Memkind memkind) {
  if (memkind == Memkind::DeviceGlobalMemory) {
#ifdef ACL_DEVICE
    auto defaultStream = device::DeviceInstance::getInstance().api->getDefaultStream();
    device::DeviceInstance::getInstance().algorithms.fillArray(
        reinterpret_cast<char*>(dst), static_cast<char>(0), size, defaultStream);
    device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
#else
    assert(false);
#endif
  } else {
    std::memset(dst, 0, size);
  }
}

void* hostToDevicePointer(void* host, enum Memkind memkind) {
  if (host == nullptr) {
    return nullptr;
  }
  if (memkind == Memkind::PinnedMemory) {
#ifdef ACL_DEVICE
    return device::DeviceInstance::getInstance().api->devicePointer(host);
#else
    return host;
#endif
  } else {
    return host;
  }
}

void printMemoryAlignment(const std::vector<std::vector<unsigned long long>>& memoryAlignment) {
  logDebug() << "Printing memory alignment per struct";
  for (unsigned long long i = 0; i < memoryAlignment.size(); i++) {
    logDebug() << memoryAlignment[i][0] << "," << memoryAlignment[i][1];
  }
}

ManagedAllocator::~ManagedAllocator() {
  for (const auto& [memkind, pointer] : dataMemoryAddresses) {
    free(pointer, memkind);
  }

  // reset memory vectors
  dataMemoryAddresses.clear();
}

void* ManagedAllocator::allocateMemory(size_t size, size_t alignment, enum Memkind memkind) {
  // NOLINTNEXTLINE(misc-const-correctness)
  void* const ptrBuffer = allocate(size, alignment, memkind);
  dataMemoryAddresses.emplace_back(memkind, ptrBuffer);
  return ptrBuffer;
}

} // namespace seissol::memory
