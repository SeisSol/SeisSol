// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_COMMON_SEGMENTMAP_H_
#define SEISSOL_SRC_COMMON_SEGMENTMAP_H_

#include <cmath>
#include <map>
#include <numeric>
#include <optional>
#include <utility>
#include <utils/logger.h>

namespace seissol {

/**
    A map that supports range inputs and point queries.
    Somewhat remniscient to a segment tree, but without a fancy query function.

    How it works internally: use the tree structure of a std::map.
    Add the start and end range points for each value.
    Add sentinel (i.e. empty optional) ranges to fill the rest.

    Query via upper_bound and lower_bound
 */
template <typename KeyT, typename ValueT>
class SegmentMap {
  public:
  SegmentMap() {
    ranges_[std::numeric_limits<KeyT>::min()] = {};
    ranges_[std::numeric_limits<KeyT>::max()] = {};
  }

  /**
      Inserts a range. start and end are both inclusive; but can be unbounded (== empty optional).
   */
  void addRange(std::optional<KeyT> start, std::optional<KeyT> end, ValueT type) {
    const auto trueStart = start.value_or(std::numeric_limits<KeyT>::min());
    const auto trueEnd = end.value_or(std::numeric_limits<KeyT>::max());

    if (trueStart > trueEnd) {
      logError() << "Invalid range:" << trueStart << "to" << trueEnd;
    }

    // check if we're empty first (i.e. look for a lower bound; then increment once; we should be
    // outside the interval again)
    const auto lb = ranges_.lower_bound(trueEnd);
    const auto pb = std::prev(lb);

    if ((lb->first == trueEnd && lb->second.has_value()) || pb->first > trueStart ||
        (pb->first == trueStart && pb->second.has_value())) {
      logError() << "Map error: range not empty." << trueStart << "to" << trueEnd;
      throw;
    }

    // insert the ranges and sentinel ranges. Also take care of any overflows/underflows that could
    // occur. (in case we would not have a sentinel range)

    ranges_[trueEnd] = type;
    if (trueEnd < std::numeric_limits<KeyT>::max()) {
      const auto pos = trueEnd + 1;
      ranges_.try_emplace(pos);
    }

    ranges_[trueStart] = type;
    if (trueStart > std::numeric_limits<KeyT>::min()) {
      const auto pos = trueStart - 1;
      ranges_.try_emplace(pos);
    }
  }

  /**
    Add a point entry.
    */
  void addEntry(KeyT value, ValueT type) { addRange(value, value, type); }

  /**
    Query a value in the segment map. Will return an empty optional if nothing was found.
   */
  [[nodiscard]] std::optional<ValueT> at(KeyT index) const {
    // (ab)use a std::map as a sort of segment tree; to support ranges
    const auto exact = ranges_.find(index);
    const auto upper = ranges_.upper_bound(index);
    const auto lower = ranges_.lower_bound(index);

    if (exact != ranges_.end()) {
      return exact->second;
    } else if (upper->second == lower->second) {
      return lower->second;
    } else {
      // error (we should not end up here)
      logError() << "Map error: internal issue.";
      throw;
    }
  }

  private:
  std::map<KeyT, std::optional<ValueT>> ranges_;
};

} // namespace seissol
#endif // SEISSOL_SRC_COMMON_SEGMENTMAP_H_
