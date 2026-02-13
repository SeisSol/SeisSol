// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_COMPACTOPTIONAL_H_
#define SEISSOL_SRC_COMMON_COMPACTOPTIONAL_H_

#include <limits>
#include <optional>
#include <stdexcept>

namespace seissol {

/**

  An optional implementation that does not consume extra memory.
  However, you'll need to specify one empty value.

 */
template <typename T, T Empty>
class CompactOptional {
  private:
  T value_{Empty};

  public:
  // NOLINTNEXTLINE
  CompactOptional(T value) : value_(value) {}
  CompactOptional() = default;

  [[nodiscard]] bool hasValue() const { return value_ != Empty; }

  [[nodiscard]] T value() const {
    if (!hasValue()) {
      throw std::runtime_error("The optional has no value, but it was requested.");
    }
    return value_;
  }

  [[nodiscard]] T valueOr(T alternative) const {
    if (hasValue()) {
      return value_;
    } else {
      return alternative;
    }
  }

  [[nodiscard]] std::optional<T> asOptional() const {
    if (hasValue()) {
      return std::optional<T>{value_};
    } else {
      return std::optional<T>{};
    }
  }
};

using OptionalSize = CompactOptional<std::size_t, std::numeric_limits<std::size_t>::max()>;

} // namespace seissol

#endif // SEISSOL_SRC_COMMON_COMPACTOPTIONAL_H_
