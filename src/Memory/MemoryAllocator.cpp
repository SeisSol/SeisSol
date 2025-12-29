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

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <stdlib.h>
#include <utils/logger.h>
#include <vector>

#ifdef ACL_DEVICE
#include "device.h"
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
  std::size_t trueSize = ((size + alignment - 1) / alignment) * alignment;
  std::size_t totalSize = trueSize;

  std::size_t start = 0;

  MPI_Allreduce(&trueSize,
                &totalSize,
                1,
                seissol::MPI::castToMpiType<std::size_t>(),
                MPI_SUM,
                seissol::MPI::mpi.comm());
  MPI_Exscan(&trueSize,
             &start,
             1,
             seissol::MPI::castToMpiType<std::size_t>(),
             MPI_SUM,
             seissol::MPI::mpi.comm());

  void* ptr = mallocShmemInternal(totalSize, alignment);

  uint8_t* dataPtr = reinterpret_cast<uint8_t*>(ptr);
  dataPtr += start;

  return reinterpret_cast<void*>(dataPtr);
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

void freeShmem(void* ptr) {
  void* startPtr = ptr;
  MPI_Bcast(&startPtr, 1, seissol::MPI::castToMpiType<std::size_t>(), 0, seissol::MPI::mpi.comm());

  freeShmemInternal(startPtr);
}

} // namespace

#endif

namespace seissol::memory {

void* allocate(size_t size, size_t alignment, enum Memkind memkind) {
  // presumably a clang-tidy bug; it'd want to make ptrBuffer const
  // NOLINTNEXTLINE(misc-const-correctness)
  void* ptrBuffer{nullptr};
  bool error = false;

  /* handle zero allocation */
  if (size == 0 && memkind != Memkind::Shmem) {
    // logWarning() << "allocation of size 0 requested, returning nullptr; (alignment: " <<
    // alignment << ", memkind: " << memkind << ").";
    ptrBuffer = nullptr;
    return ptrBuffer;
  }

#if defined(USE_MEMKIND) || defined(ACL_DEVICE)
  if (memkind == 0) {
#endif
    if (alignment % (sizeof(void*)) != 0) {
      ptrBuffer = malloc(size);
      error = (ptrBuffer == nullptr);
    } else {
      error = (posix_memalign(&ptrBuffer, alignment, size) != 0);
    }
#ifdef USE_MEMKIND
  } else if (memkind == HighBandwidth) {
    if (alignment % (sizeof(void*)) != 0) {
      ptrBuffer = hbw_malloc(size);
      error = (ptrBuffer == nullptr);
    } else {
      error = (hbw_posix_memalign(&ptrBuffer, alignment, size) != 0);
    }
#endif

#ifdef ACL_DEVICE
  } else if (memkind == DeviceGlobalMemory) {
    ptrBuffer = device::DeviceInstance::getInstance().api->allocGlobMem(size);
  } else if (memkind == DeviceUnifiedMemory) {
    ptrBuffer = device::DeviceInstance::getInstance().api->allocUnifiedMem(size);
  } else if (memkind == PinnedMemory) {
    ptrBuffer = device::DeviceInstance::getInstance().api->allocPinnedMem(size);
  } else if (memkind == DeviceGlobalCompressed) {
    ptrBuffer = device::DeviceInstance::getInstance().api->allocGlobMem(size, true);
#endif

#ifdef USE_SHMEM
  } else if (memkind == Memkind::Shmem) {
    ptrBuffer = mallocShmem(size);
#endif

#if defined(USE_MEMKIND) || defined(ACL_DEVICE)
  } else {
    logError() << "unknown memkind type used (" << memkind
               << "). Please, refer to the documentation";
  }
#endif

  if (error) {
    logError() << "The malloc failed (bytes: " << size << ", alignment: " << alignment
               << ", memkind: " << memkind << ").";
  }

  return ptrBuffer;
}

void free(void* pointer, enum Memkind memkind) {
#if defined(USE_MEMKIND) || defined(ACL_DEVICE)
  if (memkind == Standard) {
#endif
    ::free(pointer);
#ifdef USE_MEMKIND
  } else if (memkind == HighBandwidth) {
    hbw_free(pointer);
#endif

#ifdef ACL_DEVICE
  } else if (memkind == DeviceGlobalMemory) {
    device::DeviceInstance::getInstance().api->freeGlobMem(pointer);
  } else if (memkind == DeviceUnifiedMemory) {
    device::DeviceInstance::getInstance().api->freeUnifiedMem(pointer);
  } else if (memkind == PinnedMemory) {
    device::DeviceInstance::getInstance().api->freePinnedMem(pointer);
  } else if (memkind == DeviceGlobalCompressed) {
    device::DeviceInstance::getInstance().api->freeGlobMem(pointer);
#endif

#ifdef USE_SHMEM
  } else if (memkind == Memkind::Shmem) {
    freeShmem(pointer);
#endif

#if defined(USE_MEMKIND) || defined(ACL_DEVICE)
  } else {
    logError() << "unknown memkind type used (" << memkind
               << "). Please, refer to the documentation";
  }
#endif
}

void memcopy(void* dst,
             const void* src,
             std::size_t size,
             enum Memkind dstMemkind,
             enum Memkind srcMemkind) {
  if (dstMemkind == DeviceGlobalMemory && srcMemkind != DeviceGlobalMemory) {
#ifdef ACL_DEVICE
    device::DeviceInstance::getInstance().api->copyTo(dst, src, size);
#endif
  } else if (dstMemkind != DeviceGlobalMemory && srcMemkind == DeviceGlobalMemory) {
#ifdef ACL_DEVICE
    device::DeviceInstance::getInstance().api->copyFrom(dst, src, size);
#endif
  } else if (dstMemkind == DeviceGlobalMemory && srcMemkind == DeviceGlobalMemory) {
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

void printMemoryAlignment(std::vector<std::vector<unsigned long long>> memoryAlignment) {
  logDebug() << "printing memory alignment per struct";
  for (unsigned long long i = 0; i < memoryAlignment.size(); i++) {
    logDebug() << memoryAlignment[i][0] << ", " << memoryAlignment[i][1];
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
