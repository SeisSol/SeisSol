#pragma once

#ifdef ACL_DEVICE

#include "Initializer/BasicTypedefs.hpp"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.hpp"
#include <cstddef>

namespace seissol::initializers {
struct GraphKey {
  GraphKey(ComputeGraphType userGraphType,
           double userTimeWidth = 0.0,
           bool withDisplacements = false)
      : graphType(userGraphType),
        timeWidth(userTimeWidth),
        withDisplacements(withDisplacements) {}

  bool operator==(const GraphKey& other) const {
    return ((graphType == other.graphType) &&
            (timeWidth == other.timeWidth) &&
            (withDisplacements == other.withDisplacements));
  }

  ComputeGraphType graphType{};
  double timeWidth{};
  bool withDisplacements{};
};

struct GraphKeyHash {
  std::size_t operator()(GraphKey const& key) const {
    std::size_t result = 0;
    recording::hashCombine(result, key.graphType);
    recording::hashCombine(result, key.timeWidth);
    return result;
  }
};
} // namespace seissol::initializers
#endif  // ACL_DEVICE