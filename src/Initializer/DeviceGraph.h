#pragma once

#ifdef ACL_DEVICE

#include "Initializer/BasicTypedefs.hpp"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.hpp"
#include <cstddef>

namespace seissol::initializers {
struct GraphKey {
  GraphKey(ComputeGraphType userGraphType,
           double userTimeWidth = 0.0,
           bool withDisplacements = false,
           double userSubTime = 0.0)
      : graphType(userGraphType),
        timeWidth(userTimeWidth),
        withDisplacements(withDisplacements),
        subTime(userSubTime) {}

  bool operator==(const GraphKey& other) const {
    return ((graphType == other.graphType) &&
            (timeWidth == other.timeWidth) &&
            (withDisplacements == other.withDisplacements) &&
            (subTime == other.subTime));
  }

  ComputeGraphType graphType{};
  double timeWidth{};
  bool withDisplacements{};
  double subTime{};
};

struct GraphKeyHash {
  std::size_t operator()(GraphKey const& key) const {
    std::size_t result = 0;
    recording::hashCombine(result, static_cast<size_t>(key.graphType));
    recording::hashCombine(result, key.timeWidth);
    recording::hashCombine(result, key.withDisplacements ? 1 : 0);
    recording::hashCombine(result, key.subTime);
    return result;
  }
};
} // namespace seissol::initializers
#endif  // ACL_DEVICE