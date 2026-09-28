// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_COMPACTOPTIONAL_H_
#define SEISSOL_SRC_COMMON_COMPACTOPTIONAL_H_

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <type_traits>
#include <utils/logger.h>

namespace seissol {

/**

  An optional implementation that does not consume extra memory.
  However, you'll need to specify one sentinel value.

 */
template <typename T, T Empty>
class CompactOptional {
  static_assert(std::is_trivially_copyable_v<T>,
                "CompactOptional stores its value in place; use std::optional otherwise.");

  private:
  T value_{Empty};

  public:
  // NOLINTNEXTLINE
  constexpr CompactOptional(T value) : value_(value) {}
  constexpr CompactOptional() = default;

  [[nodiscard]] constexpr bool hasValue() const noexcept { return value_ != Empty; }

  [[nodiscard]] constexpr T value() const {
    if (!hasValue()) {
      logError() << "The optional has no value, but it was requested.";
    }
    return value_;
  }

  [[nodiscard]] constexpr T valueOr(T alternative) const noexcept {
    if (hasValue()) {
      return value_;
    } else {
      return alternative;
    }
  }

  [[nodiscard]] constexpr std::optional<T> asOptional() const {
    if (hasValue()) {
      return std::optional<T>{value_};
    } else {
      return std::optional<T>{};
    }
  }

  constexpr void reset() noexcept { value_ = Empty; }

  friend constexpr bool operator==(const CompactOptional& lhs, const CompactOptional& rhs) noexcept {
    return lhs.value_ == rhs.value_;
  }

  friend constexpr bool operator!=(const CompactOptional& lhs, const CompactOptional& rhs) noexcept {
    return !(lhs == rhs);
  }
};

using OptionalSize = CompactOptional<std::size_t, std::numeric_limits<std::size_t>::max()>;

/// A side of a cell, or nothing; sides are stored as small signed integers throughout.
using OptionalSide = CompactOptional<std::int8_t, -1>;

} // namespace seissol

#endif // SEISSOL_SRC_COMMON_COMPACTOPTIONAL_H_
