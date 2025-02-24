// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_INITIALIZER_MEMORYALLOCATOR_H_
#define SEISSOL_SRC_INITIALIZER_MEMORYALLOCATOR_H_

#include "Common/Constants.h"
#include "Common/Marker.h"
#include "Kernels/Precision.h"
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
  SEISSOL_HOSTDEVICE T* begin() noexcept { return data_; }
  SEISSOL_HOSTDEVICE T* end() noexcept { return data_ + N; }
  SEISSOL_HOSTDEVICE const T* begin() const noexcept { return data_; }
  SEISSOL_HOSTDEVICE const T* end() const noexcept { return data_ + N; }
  SEISSOL_HOSTDEVICE T* data() { return data_; }
  SEISSOL_HOSTDEVICE const T* data() const { return data_; }
  SEISSOL_HOSTDEVICE constexpr T& operator[](std::size_t pos) { return data_[pos]; }
  SEISSOL_HOSTDEVICE constexpr const T& operator[](std::size_t pos) const { return data_[pos]; }
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr std::size_t size() const noexcept { return N; }

  private:
  alignas(Alignment) T data_[N];
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

void* hostToDevicePointer(void* host, enum Memkind memkind);

template <typename T>
T* hostToDevicePointerTyped(T* host, enum Memkind memkind) {
  return hostToDevicePointer(host, memkind);
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

  ManagedAllocator(const ManagedAllocator&) = delete;
  auto operator=(const ManagedAllocator&) = delete;

  ManagedAllocator(ManagedAllocator&&) = default;
  auto operator=(ManagedAllocator&&) noexcept -> ManagedAllocator& = default;

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
  MemkindArray(MemkindArray<T>&& source) = default;
  MemkindArray(const MemkindArray<T>& source) : MemkindArray(source, source.memkind) {}

  auto operator=(const MemkindArray<T>& source) -> MemkindArray& {
    if (capacity != source.capacity) {
      resize(source.capacity);
    }
    copyFrom(source);
    return *this;
  }
  auto operator=(MemkindArray<T>&& source) noexcept -> MemkindArray& = default;

  MemkindArray(const MemkindArray<T>& source, Memkind memkind)
      : MemkindArray(source.size(), memkind) {
    copyFrom(source);
  }
  MemkindArray(const std::vector<T>& source, Memkind memkind)
      : MemkindArray(source.size(), memkind) {
    copyFrom(source);
  }
  MemkindArray(std::size_t capacity, Memkind memkind) : MemkindArray(memkind) { resize(capacity); }
  MemkindArray(Memkind memkind) : memkind(memkind) {}

  void resize(std::size_t capacity) {
    this->capacity = capacity;
    free(dataPtr, memkind);
    dataPtr = allocTyped<T>(capacity, Alignment, memkind);
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
  SEISSOL_HOSTDEVICE T* data() noexcept { return dataPtr; }
  SEISSOL_HOSTDEVICE const T* data() const noexcept { return dataPtr; }
  SEISSOL_HOSTDEVICE T* begin() noexcept { return dataPtr; }
  SEISSOL_HOSTDEVICE T* end() noexcept { return dataPtr + capacity; }
  SEISSOL_HOSTDEVICE const T* begin() const noexcept { return dataPtr; }
  SEISSOL_HOSTDEVICE const T* end() const noexcept { return dataPtr + capacity; }
  SEISSOL_HOSTDEVICE constexpr T& operator[](std::size_t index) { return dataPtr[index]; }
  SEISSOL_HOSTDEVICE constexpr const T& operator[](std::size_t index) const {
    return dataPtr[index];
  }
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr std::size_t size() const noexcept { return capacity; }

  private:
  T* dataPtr{nullptr};
  std::size_t capacity{0};
  enum Memkind memkind;
};

} // namespace seissol::memory

#endif // SEISSOL_SRC_INITIALIZER_MEMORYALLOCATOR_H_
