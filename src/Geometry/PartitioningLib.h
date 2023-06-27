// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PARTITIONINGLIB_H
#define PARTITIONINGLIB_H

#include "PUML/Partition.h"
#include <string_view>

namespace seissol {

PUML::PartitionerType toPartitionerType(std::string_view partitioningLib);
std::string_view toStringView(PUML::PartitionerType type);

} // namespace seissol

#endif // PARTITIONINGLIB_H
