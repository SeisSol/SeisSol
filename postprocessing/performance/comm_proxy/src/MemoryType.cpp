// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#include "MemoryType.hpp"

namespace seissol {

std::string to_string(MemoryType type) {
  switch (type) {
  case MemoryType::Host:
    return "host";
  case MemoryType::Shared:
    return "shared";
  case MemoryType::Device:
    return "device";
  default:
    break;
  }
  return {};
}

} // namespace seissol
