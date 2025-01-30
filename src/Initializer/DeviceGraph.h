// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_DEVICEGRAPH_H_
#define SEISSOL_SRC_INITIALIZER_DEVICEGRAPH_H_

#ifdef ACL_DEVICE

#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include <cstddef>

namespace seissol::initializer {
struct GraphKey {
  GraphKey(ComputeGraphType userGraphType,
           double userTimeWidth = 0.0,
           bool withDisplacements = false)
      : graphType(userGraphType), timeWidth(userTimeWidth), withDisplacements(withDisplacements) {}

  bool operator==(const GraphKey& other) const {
    return ((graphType == other.graphType) && (timeWidth == other.timeWidth) &&
            (withDisplacements == other.withDisplacements));
  }

  ComputeGraphType graphType{};
  double timeWidth{};
  bool withDisplacements{};
};

struct GraphKeyHash {
  std::size_t operator()(const GraphKey& key) const {
    std::size_t result = 0;
    recording::hashCombine(result, static_cast<size_t>(key.graphType));
    recording::hashCombine(result, key.timeWidth);
    return result;
  }
};
} // namespace seissol::initializer

#endif // ACL_DEVICE

#endif // SEISSOL_SRC_INITIALIZER_DEVICEGRAPH_H_
