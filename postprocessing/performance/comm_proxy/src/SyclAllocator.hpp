// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SYCLALLOCATOR_20231020_HPP
#define SYCLALLOCATOR_20231020_HPP

#include "Allocator.hpp"
#include "MemoryType.hpp"

#include <sycl/sycl.hpp>
#include <cstddef>
#include <utility>

namespace seissol {

class SyclAllocator : public Allocator {
  public:
  SyclAllocator(MemoryType type, sycl::queue queue)
      : type_(usm::alloc::host), queue_(std::move(queue)) {
    switch (type) {
    case MemoryType::Host:
      type_ = usm::alloc::host;
      break;
    case MemoryType::Shared:
      type_ = usm::alloc::shared;
      break;
    case MemoryType::Device:
      type_ = usm::alloc::device;
      break;
    }
  }
  inline void* malloc(std::size_t numBytes) override { return ::malloc(numBytes); }
  inline void free(void* ptr) override { return ::free(ptr); }

  private:
  sycl::usm::alloc type_;
  sycl::queue queue_;
};

} // namespace seissol

#endif // SYCLALLOCATOR_20231020_HPP
