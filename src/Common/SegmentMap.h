// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_COMMON_SEGMENTMAP_H_
#define SEISSOL_SRC_COMMON_SEGMENTMAP_H_

#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <utils/logger.h>

namespace seissol {

/**
    A map that supports range inputs and point queries.
    Somewhat remniscient to a segment tree, but without a fancy query function.

    How it works internally: use the tree structure of a std::map.
    Add the start and end range points for each value.
    Add sentinel (i.e. empty optional) ranges to fill the rest.

    Invariants: numeric_limits<KeyT>::min() and numeric_limits<KeyT>::max() are always keys, and a
    key covers the half-open interval (predecessor key, key]. Hence lower_bound() never returns
    end() for any key, and it is the only lookup a query needs.
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

    // the first key at or above trueEnd covers trueEnd itself; a value there means we overlap
    const auto lb = ranges_.lower_bound(trueEnd);
    const auto occupiedAbove = lb->second.has_value();

    // any key inside [trueStart, trueEnd) delimits a neighboring range; it may only sit exactly on
    // trueStart, and only as a sentinel. lb == begin() means trueEnd, and hence trueStart as well,
    // is the smallest representable key -- there is nothing below it to look at.
    const auto occupiedBelow = [&]() {
      if (lb == ranges_.begin()) {
        return false;
      }
      const auto pb = std::prev(lb);
      return pb->first > trueStart || (pb->first == trueStart && pb->second.has_value());
    }();

    if (occupiedAbove || occupiedBelow) {
      logError() << "Map error: range not empty." << trueStart << "to" << trueEnd;
    }

    // insert the ranges and sentinel ranges. Also take care of any overflows/underflows that could
    // occur. (in case we would not have a sentinel range)

    ranges_[trueEnd] = type;
    if (trueEnd < std::numeric_limits<KeyT>::max()) {
      const auto pos = static_cast<KeyT>(trueEnd + 1);
      ranges_.try_emplace(pos);
    }

    ranges_[trueStart] = type;
    if (trueStart > std::numeric_limits<KeyT>::min()) {
      const auto pos = static_cast<KeyT>(trueStart - 1);
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
    return ranges_.lower_bound(index)->second;
  }

  private:
  std::map<KeyT, std::optional<ValueT>> ranges_;
};

} // namespace seissol
#endif // SEISSOL_SRC_COMMON_SEGMENTMAP_H_
