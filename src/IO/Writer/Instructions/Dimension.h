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
#include <string>
#include <utils/logger.h>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer {

//! @brief How a dimension of a data source is laid out across the ranks.
enum class Extent : std::uint8_t {
  //! Every rank holds the same entries, and they are written once.
  Replicated,
  //! The entries are split across the ranks and concatenated in rank order.
  Distributed
};

/**
 * @brief Whether a dimension is settled when the dataset is created, or grows with every write.
 *
 * Independent of how the ranks share it: a VTKHDF time series concatenates its steps along the
 * same dimension the ranks are split along, while a table of samples per receiver grows along one
 * dimension and is distributed along another.
 */
enum class Growth : std::uint8_t {
  //! The dataset has this extent once and for all.
  Fixed,
  //! The dataset is unlimited here, and every write adds its entries at the end.
  Appended
};

/**
 * @brief One dimension of a data source.
 *
 * A distributed dimension has no size here: how many entries a rank contributes is only known
 * once the data is there, and the total follows from a scan over the ranks. Every other dimension
 * carries the extent of one write -- for an appended dimension that is how far the dataset grows,
 * not how large it ends up.
 */
struct Dimension {
  std::size_t size{0};
  Extent extent{Extent::Replicated};
  Growth growth{Growth::Fixed};
  /**
   * @brief How far a storage chunk reaches along this dimension, or zero to let the writer pick.
   *
   * The chunking is settled when the dataset is created and holds for the whole run, so a
   * dimension whose per-write extent varies cannot have the writer derive it from the first
   * write and expect it to fit the rest. Which value pays off depends on how the file is read
   * afterwards, which is not something the writer can know.
   */
  std::size_t chunk{0};

  static Dimension replicated(std::size_t size) {
    return {size, Extent::Replicated, Growth::Fixed, 0};
  }
  static Dimension distributed() { return {0, Extent::Distributed, Growth::Fixed, 0}; }

  /**
   * @brief A dimension every write extends by @p size entries, the same number on every rank.
   *
   * @p chunk is the storage chunk extent along it; zero takes @p size, which is right where
   * every write contributes the same amount.
   */
  static Dimension appended(std::size_t size, std::size_t chunk = 0) {
    return {size, Extent::Replicated, Growth::Appended, chunk};
  }

  /**
   * @brief A dimension that is both split across the ranks and extended by every write.
   *
   * What a flat array a reader slices by itself looks like: one write contributes a round of
   * distributed data, and the next one lands behind it.
   */
  static Dimension distributedAppended() { return {0, Extent::Distributed, Growth::Appended, 0}; }

  [[nodiscard]] bool isDistributed() const { return extent == Extent::Distributed; }
  [[nodiscard]] bool isAppended() const { return growth == Growth::Appended; }
};

/**
 * @brief Builds a shape from fixed sizes, optionally preceded by a distributed dimension.
 *
 * The common case, and the one the mesh writers use: one dimension per rank-local entry, and
 * fixed trailing dimensions describing what each entry holds.
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

/**
 * @brief Rejects shapes the writers cannot express, naming what is wrong.
 *
 * A second distributed dimension would make the offset of a rank depend on more than one scan,
 * and a second appended one would leave open which of them a write extends.
 */
inline void checkDimensions(const std::vector<Dimension>& dimensions, const std::string& what) {
  std::size_t distributed = 0;
  std::size_t appended = 0;
  for (const auto& dimension : dimensions) {
    distributed += static_cast<std::size_t>(dimension.isDistributed());
    appended += static_cast<std::size_t>(dimension.isAppended());
  }
  if (distributed > 1) {
    logError() << "The shape of" << what << "has" << distributed
               << "distributed dimensions; at most one is supported.";
  }
  if (appended > 1) {
    logError() << "The shape of" << what << "has" << appended
               << "appended dimensions; at most one is supported.";
  }
}

inline YAML::Node serializeDimensions(const std::vector<Dimension>& dimensions) {
  YAML::Node node;
  for (const auto& dimension : dimensions) {
    YAML::Node entry;
    entry["size"] = dimension.size;
    entry["distributed"] = dimension.isDistributed();
    entry["appended"] = dimension.isAppended();
    entry["chunk"] = dimension.chunk;
    node.push_back(entry);
  }
  return node;
}

inline std::vector<Dimension> deserializeDimensions(const YAML::Node& node) {
  std::vector<Dimension> result;
  result.reserve(node.size());
  for (const auto& entry : node) {
    result.push_back(
        Dimension{entry["size"].as<std::size_t>(),
                  entry["distributed"].as<bool>() ? Extent::Distributed : Extent::Replicated,
                  entry["appended"].as<bool>() ? Growth::Appended : Growth::Fixed,
                  entry["chunk"].as<std::size_t>()});
  }
  return result;
}

} // namespace seissol::io::writer

#endif // SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_DIMENSION_H_
