// SPDX-FileCopyrightText: 2024 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_GEOMETRY_PARTITIONINGLIB_H_
#define SEISSOL_SRC_GEOMETRY_PARTITIONINGLIB_H_

#include "PUML/Partition.h"
#include <string_view>

namespace seissol {

PUML::PartitionerType toPartitionerType(std::string_view partitioningLib);
std::string_view toStringView(PUML::PartitionerType type);

} // namespace seissol

#endif // SEISSOL_SRC_GEOMETRY_PARTITIONINGLIB_H_
