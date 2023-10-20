// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#include <string>

namespace seissol {

enum class MemoryType { Host, Shared, Device };

std::string to_string(MemoryType type);

} // namespace seissol
