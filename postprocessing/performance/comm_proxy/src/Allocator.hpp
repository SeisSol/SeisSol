// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef ALLOCATOR_20231020_HPP
#define ALLOCATOR_20231020_HPP

#include <cstddef>

namespace seissol {

class Allocator {
  public:
  virtual ~Allocator() {}

  virtual void* malloc(std::size_t numBytes) = 0;
  virtual void free(void* ptr) = 0;
};

} // namespace seissol

#endif // ALLOCATOR_20231020_HPP
