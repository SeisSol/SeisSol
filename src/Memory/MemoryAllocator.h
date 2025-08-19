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

#ifndef SEISSOL_SRC_MEMORY_MEMORYALLOCATOR_H_
#define SEISSOL_SRC_MEMORY_MEMORYALLOCATOR_H_

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
#include <Device/UsmAllocator.h>
#include <Device/device.h>
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
  [[nodiscard]] SEISSOL_HOSTDEVICE const T* begin() const noexcept { return data_; }
  [[nodiscard]] SEISSOL_HOSTDEVICE const T* end() const noexcept { return data_ + N; }
  SEISSOL_HOSTDEVICE T* data() { return data_; }
  [[nodiscard]] SEISSOL_HOSTDEVICE const T* data() const { return data_; }
  SEISSOL_HOSTDEVICE constexpr T& operator[](std::size_t pos) { return data_[pos]; }
  SEISSOL_HOSTDEVICE constexpr const T& operator[](std::size_t pos) const { return data_[pos]; }
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr std::size_t size() const noexcept { return N; }

  private:
  alignas(Alignment) T data_[N];
};

enum class Memkind {
  Standard = 0,
  HighBandwidth = 1,
  DeviceGlobalMemory = 3,
  DeviceUnifiedMemory = 4,
  PinnedMemory = 5
};

void* allocate(size_t size, size_t alignment = 1, Memkind memkind = Memkind::Standard);

template <typename T>
T* allocTyped(size_t count, size_t alignment = 1, Memkind memkind = Memkind::Standard) {
  return reinterpret_cast<T*>(allocate(count * sizeof(T), alignment, memkind));
}

void free(void* pointer, Memkind memkind = Memkind::Standard);

void memcopy(
    void* dst, const void* src, std::size_t size, enum Memkind dstMemkind, enum Memkind srcMemkind);

template <typename T>
void memcopyTyped(
    T* dst, const T* src, std::size_t count, enum Memkind dstMemkind, enum Memkind srcMemkind) {
  // explicit cast to silence clang-tidy
  memcopy(static_cast<void*>(dst),
          static_cast<const void*>(src),
          count * sizeof(T),
          dstMemkind,
          srcMemkind);
}

void memzero(void* dst, std::size_t size, enum Memkind memkind);

template <typename T>
void memzeroTyped(T* dst, std::size_t count, enum Memkind memkind) {
  memzero(dst, count * sizeof(T), memkind);
}

template <typename T>
void meminit(T* dst, std::size_t count, enum Memkind memkind) {
  std::vector<T> local(count);
  memcopyTyped<T>(dst, local.data(), count, memkind, Memkind::Standard);
}

void* hostToDevicePointer(void* host, enum Memkind memkind);

template <typename T>
T* hostToDevicePointerTyped(T* host, enum Memkind memkind) {
  return reinterpret_cast<T*>(hostToDevicePointer(host, memkind));
}

/**
 * Prints the memory alignment of in terms of relative start and ends in bytes.
 *
 * @param memoryAlignment memory alignment.
 **/
void printMemoryAlignment(const std::vector<std::vector<unsigned long long>>& memoryAlignment);

/**
 * Automatically frees allocated memory on destruction.
 **/
class ManagedAllocator {
  private:
  using Address = std::pair<enum Memkind, void*>;
  using AddressVector = std::vector<Address>;

  //! holds all memory addresses, which point to data arrays and have been returned by mallocs
  //! calling functions of the memory allocator.
  AddressVector dataMemoryAddresses_;

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
  void* allocateMemory(size_t size, size_t alignment = 1, Memkind memkind = Memkind::Standard);
};

template <typename T>
class MemkindArray {
  public:
  MemkindArray(MemkindArray<T>&& source) noexcept
      : dataPtr_(source.dataPtr_), capacity_(source.capacity_), memkind_(source.memkind_) {
    source.dataPtr_ = nullptr;
    source.memkind_ = Memkind::Standard;
    source.capacity_ = 0;
  }
  MemkindArray(const MemkindArray<T>& source) : MemkindArray(source, source.memkind_) {}

  auto operator=(const MemkindArray<T>& source) -> MemkindArray& {
    if (&source != this) {
      if (capacity_ != source.capacity_ || memkind_ != source.memkind_) {
        resize(source.capacity_, source.memkind_);
      }
      copyFrom(source);
    }
    return *this;
  }
  auto operator=(MemkindArray<T>&& source) noexcept -> MemkindArray& {
    if (&source != this) {
      this->dataPtr_ = source.dataPtr_;
      this->memkind_ = source.memkind_;
      this->capacity_ = source.capacity_;

      source.dataPtr_ = nullptr;
      source.memkind_ = Memkind::Standard;
      source.capacity_ = 0;
    }
    return *this;
  }

  MemkindArray(const MemkindArray<T>& source, Memkind memkind)
      : MemkindArray(source.size(), memkind) {
    copyFrom(source);
  }
  MemkindArray(const std::vector<T>& source, Memkind memkind)
      : MemkindArray(source.size(), memkind) {
    copyFrom(source);
  }
  MemkindArray(std::size_t capacity, Memkind memkind) : MemkindArray(memkind) {
    resize(capacity, memkind);
  }
  explicit MemkindArray(Memkind memkind) : memkind_(memkind) {}

  void resize(std::size_t capacity, Memkind newMemkind) {
    this->capacity_ = capacity;
    free(dataPtr_, memkind_);
    dataPtr_ = allocTyped<T>(capacity, Alignment, newMemkind);
    this->memkind_ = newMemkind;
  }
  void resize(std::size_t capacity) { resize(capacity, memkind_); }

  void copyFrom(const std::vector<T>& source) {
    assert(source.size() <= capacity);
    memcopyTyped<T>(dataPtr_, source.data(), capacity_, memkind_, Memkind::Standard);
  }
  void copyFrom(const MemkindArray<T>& source) {
    assert(source.size() <= capacity);
    memcopyTyped<T>(dataPtr_, source.data(), capacity_, memkind_, source.memkind_);
  }
  ~MemkindArray() { free(dataPtr_, memkind_); }
  SEISSOL_HOSTDEVICE T* data() noexcept { return dataPtr_; }
  [[nodiscard]] SEISSOL_HOSTDEVICE const T* data() const noexcept { return dataPtr_; }
  SEISSOL_HOSTDEVICE T* begin() noexcept { return dataPtr_; }
  SEISSOL_HOSTDEVICE T* end() noexcept { return dataPtr_ + capacity_; }
  [[nodiscard]] SEISSOL_HOSTDEVICE const T* begin() const noexcept { return dataPtr_; }
  [[nodiscard]] SEISSOL_HOSTDEVICE const T* end() const noexcept { return dataPtr_ + capacity_; }
  SEISSOL_HOSTDEVICE constexpr T& operator[](std::size_t index) { return dataPtr_[index]; }
  SEISSOL_HOSTDEVICE constexpr const T& operator[](std::size_t index) const {
    return dataPtr_[index];
  }
  [[nodiscard]] SEISSOL_HOSTDEVICE constexpr std::size_t size() const noexcept { return capacity_; }

  private:
  T* dataPtr_{nullptr};
  std::size_t capacity_{0};
  Memkind memkind_{Memkind::Standard};
};

} // namespace seissol::memory

#endif // SEISSOL_SRC_MEMORY_MEMORYALLOCATOR_H_
