// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_DIMENSION_H_
#define SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_DIMENSION_H_

#include <cstddef>
#include <cstdint>
#include <vector>

namespace seissol::io::writer {

//! @brief How a dimension of a data source is laid out across the ranks.
enum class Extent : std::uint8_t {
  //! Every rank holds the same entries, and they are written once.
  Replicated,
  //! The entries are split across the ranks and concatenated in rank order.
  Distributed
};

/**
 * @brief One dimension of a data source.
 *
 * A distributed dimension has no size here: how many entries a rank contributes is only known
 * once the data is there, and the total follows from a scan over the ranks.
 */
struct Dimension {
  std::size_t size{0};
  Extent extent{Extent::Replicated};

  static Dimension replicated(std::size_t size) { return {size, Extent::Replicated}; }
  static Dimension distributed() { return {0, Extent::Distributed}; }

  [[nodiscard]] bool isDistributed() const { return extent == Extent::Distributed; }
};

/**
 * @brief Builds a shape from replicated sizes, optionally preceded by a distributed dimension.
 *
 * At most one dimension may be distributed, and it has to be the first one. Nothing in the
 * writers needs more than that, while an arbitrary placement would turn the offset of a rank
 * from a single scan into a per-dimension one.
 */
inline std::vector<Dimension> makeDimensions(const std::vector<std::size_t>& replicated,
                                             bool leadingDistributed) {
  std::vector<Dimension> result;
  result.reserve(replicated.size() + (leadingDistributed ? 1 : 0));
  if (leadingDistributed) {
    result.push_back(Dimension::distributed());
  }
  for (const auto size : replicated) {
    result.push_back(Dimension::replicated(size));
  }
  return result;
}

} // namespace seissol::io::writer

#endif // SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_DIMENSION_H_
