// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HOSTALLOCATOR_20231020_HPP
#define HOSTALLOCATOR_20231020_HPP

#include "Allocator.hpp"

#include <cstddef>
#include <cstdio>

namespace seissol {

class HostAllocator : public Allocator {
  public:
  inline void* malloc(std::size_t numBytes) override { return ::malloc(numBytes); }
  inline void free(void* ptr) override { return ::free(ptr); }
};

} // namespace seissol

#endif // HOSTALLOCATOR_20231020_HPP
