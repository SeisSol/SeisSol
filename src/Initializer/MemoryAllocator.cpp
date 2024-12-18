/******************************************************************************
** Copyright (c) 2015, Intel Corporation                                     **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Alexander Heinecke (Intel Corp.)
******************************************************************************/
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Aligned memory allocation.
 **/
#include "MemoryAllocator.h"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <stdlib.h>
#include <utils/logger.h>
#include <omp.h>

#ifdef ACL_DEVICE
#include "device.h"
#endif

namespace seissol::memory {

void* allocate(size_t size, size_t alignment, enum Memkind memkind) {
  void* ptrBuffer{nullptr};
  bool error = false;

  // enforce USM
  memkind = DeviceUnifiedMemory;

  /* handle zero allocation */
  if (size == 0) {
    // logWarning() << "allocation of size 0 requested, returning nullptr; (alignment: " <<
    // alignment << ", memkind: " << memkind << ").";
    ptrBuffer = nullptr;
    return ptrBuffer;
  }

#if defined(USE_MEMKIND) || defined(ACL_DEVICE)
  if (memkind == 0) {
#endif
      if (alignment % (sizeof(void*)) != 0) {
        ptrBuffer = omp_alloc( size );
        error = (ptrBuffer == nullptr);
      } else {
        // error = (posix_memalign( &ptrBuffer, alignment, size ) != 0);
        ptrBuffer = omp_aligned_alloc(alignment, size);
        error = (ptrBuffer == nullptr);
      }
#ifdef USE_MEMKIND
    }
#endif

#ifdef ACL_DEVICE
  } else if (memkind == DeviceGlobalMemory) {
    ptrBuffer = device::DeviceInstance::getInstance().api->allocGlobMem(size);
  } else if (memkind == DeviceUnifiedMemory) {
    ptrBuffer = device::DeviceInstance::getInstance().api->allocUnifiedMem(size);
  } else if (memkind == PinnedMemory) {
    ptrBuffer = device::DeviceInstance::getInstance().api->allocPinnedMem(size);
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
  memkind = DeviceUnifiedMemory;
  
#if defined(USE_MEMKIND) || defined(ACL_DEVICE)
  if (memkind == Standard) {
#endif
    ::omp_free(pointer);
#ifdef USE_MEMKIND
  } else if (memkind == HighBandwidth) {
    hbw_free(pointer);
#endif

#ifdef ACL_DEVICE
  } else if ((memkind == DeviceGlobalMemory) || (memkind == DeviceUnifiedMemory)) {
    device::DeviceInstance::getInstance().api->freeMem(pointer);
  } else if (memkind == PinnedMemory) {
    device::DeviceInstance::getInstance().api->freePinnedMem(pointer);
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
  void* ptrBuffer = allocate(size, alignment, memkind);
  dataMemoryAddresses.emplace_back(memkind, ptrBuffer);
  return ptrBuffer;
}

} // namespace seissol::memory
