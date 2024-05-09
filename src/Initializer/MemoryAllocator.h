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

#ifndef MEMORYALLOCATOR_H_
#define MEMORYALLOCATOR_H_

#include "Kernels/precision.hpp"
#include <cassert>
#include <cstdlib>
#include <vector>

#ifdef USE_MEMKIND
#include <hbwmalloc.h>
#endif

#ifdef ACL_DEVICE
#include <UsmAllocator.h>
#include <device.h>
#endif

namespace seissol::memory {
#ifdef ACL_DEVICE
using AllocatorT = device::UsmAllocator<real>;
#else
using AllocatorT = std::allocator<real>;
#endif
template <typename T>
using VectorT =
    std::vector<T, typename std::allocator_traits<AllocatorT>::template rebind_alloc<T>>;

template <typename T, std::size_t N>
class AlignedArray {
  public:
  inline T* data() { return data_; }
  inline const T* data() const { return data_; }
  constexpr T& operator[](std::size_t pos) { return data_[pos]; }
  constexpr const T& operator[](std::size_t pos) const { return data_[pos]; }
  constexpr std::size_t size() const noexcept { return N; }

  private:
  alignas(ALIGNMENT) T data_[N];
};

// TODO(David): make enum class
enum Memkind {
  Standard = 0,
  HighBandwidth = 1,
  DeviceGlobalMemory = 3,
  DeviceUnifiedMemory = 4,
  PinnedMemory = 5
};
void* allocate(size_t size, size_t alignment = 1, enum Memkind memkind = Standard);

template <typename T>
T* allocTyped(size_t count, size_t alignment = 1, enum Memkind memkind = Standard) {
  return reinterpret_cast<T*>(allocate(count * sizeof(T), alignment, memkind));
}

void free(void* pointer, enum Memkind memkind = Standard);

void memcopy(
    void* dst, const void* src, std::size_t size, enum Memkind dstMemkind, enum Memkind srcMemkind);

template <typename T>
void memcopyTyped(
    T* dst, const T* src, std::size_t count, enum Memkind dstMemkind, enum Memkind srcMemkind) {
  memcopy(dst, src, count * sizeof(T), dstMemkind, srcMemkind);
}

void memzero(void* dst, std::size_t size, enum Memkind memkind);

template <typename T>
void memzeroTyped(T* dst, std::size_t count, enum Memkind memkind) {
  memzero(dst, count * sizeof(T), memkind);
}

template <typename T>
void meminit(T* dst, std::size_t count, enum Memkind memkind) {
  std::vector<T> local(count);
  memcopyTyped<T>(dst, local, count, memkind, Memkind::Standard);
}

/**
 * Prints the memory alignment of in terms of relative start and ends in bytes.
 *
 * @param memoryAlignment memory alignment.
 **/
void printMemoryAlignment(std::vector<std::vector<unsigned long long>> memoryAlignment);

/**
 * Automatically frees allocated memory on destruction.
 **/
class ManagedAllocator {
  private:
  using Address = std::pair<enum Memkind, void*>;
  using AddressVector = std::vector<Address>;

  //! holds all memory addresses, which point to data arrays and have been returned by mallocs
  //! calling functions of the memory allocator.
  AddressVector dataMemoryAddresses;

  public:
  ManagedAllocator() = default;

  /**
   * Frees all memory, which was allocated by functions of the ManagedAllocator.
   **/
  ~ManagedAllocator();

  /**
   * Allocates a single chunk of memory with the given size and alignment.
   *
   * @param  size size of the chunk in byte.
   * @param  alignment alignment of the memory chunk in byte.
   * @return pointer, which points to the aligned memory of the given size.
   **/
  void* allocateMemory(size_t size, size_t alignment = 1, enum Memkind memkind = Standard);
};

template <typename T>
class MemkindArray {
  public:
  MemkindArray(const MemkindArray<T>& source) : MemkindArray(source, source.memkind) {}
  MemkindArray(const MemkindArray<T>& source, Memkind memkind)
      : MemkindArray(source.size(), memkind) {
    copyFrom(source);
  }
  MemkindArray(const std::vector<T>& source, Memkind memkind)
      : MemkindArray(source.size(), memkind) {
    copyFrom(source);
  }
  MemkindArray(std::size_t capacity, Memkind memkind) : MemkindArray(memkind) { resize(capacity); }
  MemkindArray(Memkind memkind) : memkind(memkind), dataPtr(nullptr), capacity(0) {}
  void resize(std::size_t capacity) {
    this->capacity = capacity;
    free(dataPtr, memkind);
    dataPtr = allocTyped<T>(capacity, ALIGNMENT, memkind);
  }
  void copyFrom(const std::vector<T>& source) {
    assert(source.size() <= capacity);
    memcopyTyped<T>(dataPtr, source.data(), capacity, memkind, Memkind::Standard);
  }
  void copyFrom(const MemkindArray<T>& source) {
    assert(source.size() <= capacity);
    memcopyTyped<T>(dataPtr, source.data(), capacity, memkind, source.memkind);
  }
  ~MemkindArray() { free(dataPtr, memkind); }
  inline T* data() noexcept { return dataPtr; }
  inline const T* data() const noexcept { return dataPtr; }
  inline T* begin() noexcept { return dataPtr; }
  inline T* end() noexcept { return dataPtr + capacity; }
  inline const T* begin() const noexcept { return dataPtr; }
  inline const T* end() const noexcept { return dataPtr + capacity; }
  constexpr T& operator[](std::size_t index) { return dataPtr[index]; }
  constexpr const T& operator[](std::size_t index) const { return dataPtr[index]; }
  constexpr std::size_t size() const noexcept { return capacity; }

  private:
  T* dataPtr{nullptr};
  std::size_t capacity;
  enum Memkind memkind;
};

} // namespace seissol::memory

#endif
