// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_MEMORY_TREE_BACKMAP_H_
#define SEISSOL_SRC_MEMORY_TREE_BACKMAP_H_

#include <array>
#include <cassert>
#include <cstddef>
#include <limits>
#include <optional>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer {

/**
  A position in the Storage structure; identified by the layer color and the cell position inside
  that layer. (and the global ID as extra info)
 */
struct StoragePosition {
  std::size_t color{std::numeric_limits<std::size_t>::max()};
  std::size_t cell{std::numeric_limits<std::size_t>::max()};
  std::size_t global{std::numeric_limits<std::size_t>::max()};

  bool operator==(const StoragePosition& other) const {
    return color == other.color && cell == other.cell;
  }

  bool operator!=(const StoragePosition& other) const { return !(*this == other); }

  const static StoragePosition NullPosition;

  StoragePosition() = default;

  StoragePosition(std::size_t color, std::size_t cell, std::size_t global)
      : color(color), cell(cell), global(global) {}
};

/**
  A undefined storage position; e.g. for use when there is no cell present but a position needs to
  be stored. Will most likely cause a segfault when used; prefer using optionals to StoragePositions
  instead where possible.
 */
const inline StoragePosition StoragePosition::NullPosition =
    StoragePosition{std::numeric_limits<std::size_t>::max(),
                    std::numeric_limits<std::size_t>::max(),
                    std::numeric_limits<std::size_t>::max()};

/**
  A map from e.g. cell ID to position in a layered Storage structure. May contain MaxDuplicatesP
  many duplicates per cell ID. (akin to the Lut implementation in older SeisSol)
 */
template <std::size_t MaxDuplicatesP>
class StorageBackmap {
  public:
  static constexpr std::size_t MaxDuplicates = MaxDuplicatesP;

  private:
  using CellStoragePosition = std::array<StoragePosition, MaxDuplicates>;

  static StoragePosition getNullPosition(std::size_t index) {
    return StoragePosition::NullPosition;
  }

  template <std::size_t... Indices>
  static CellStoragePosition getNullStoragePosition(std::index_sequence<Indices...> /*...*/) {
    return CellStoragePosition{getNullPosition(Indices)...};
  }

  const inline static std::vector<StoragePosition> PreNullCellStoragePosition =
      std::vector<StoragePosition>(MaxDuplicates, StoragePosition::NullPosition);

  const inline static CellStoragePosition NullCellStoragePosition =
      getNullStoragePosition(std::make_index_sequence<MaxDuplicates>());

  std::vector<CellStoragePosition> data;

  public:
  StorageBackmap(std::size_t size) : data(NullCellStoragePosition, size) {}

  StorageBackmap() = default;

  [[nodiscard]] StoragePosition get(std::size_t id) const {
    const auto position = getDup(id, 0);
    assert(position.has_value());
    return position.value();
  }

  [[nodiscard]] std::optional<StoragePosition> getDup(std::size_t id, std::size_t duplicate) const {
    if (duplicate >= MaxDuplicates) {
      return {};
    } else {
      const auto position = data[id][duplicate];
      if (position == StoragePosition::NullPosition) {
        return {};
      } else {
        return position;
      }
    }
  }

  void setSize(std::size_t size) { data.resize(size, NullCellStoragePosition); }

  template <typename TRef>
  std::size_t addElement(std::size_t color,
                         const TRef* zeroPosition,
                         const TRef* layerData,
                         std::size_t cell,
                         std::size_t index) {
    const auto layerPosition = (reinterpret_cast<std::intptr_t>(layerData) -
                                reinterpret_cast<std::intptr_t>(zeroPosition)) /
                               sizeof(TRef);
    const auto position = layerPosition + index;
    const auto storagePosition = StoragePosition{color, index, position};
    for (std::size_t j = 0; j < MaxDuplicates; ++j) {
      if (data[cell][j] == StoragePosition::NullPosition) {
        data[cell][j] = storagePosition;
        return j;
      }
    }
    logError() << "Backmap has stored more than" << MaxDuplicates << "duplicates for" << cell
               << ". Out of capacity.";
    throw;
  }
};

} // namespace seissol::initializer
#endif // SEISSOL_SRC_MEMORY_TREE_BACKMAP_H_
