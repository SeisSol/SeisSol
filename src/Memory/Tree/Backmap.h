// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_MEMORY_TREE_BACKMAP_H_
#define SEISSOL_SRC_MEMORY_TREE_BACKMAP_H_

#include <Memory/Tree/LTSTree.h>
#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

namespace seissol::initializer {
template <std::size_t MaxDuplicates, typename LtsTree>
class StorageBackmap {
  public:
  using StoragePosition = std::pair<std::size_t, std::size_t>;
  const inline static StoragePosition NullPosition = StoragePosition(
      std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max());

  private:
  using CellStoragePosition = std::array<StoragePosition, MaxDuplicates>;

  const inline static std::vector<StoragePosition> PreNullCellStoragePosition =
      std::vector<StoragePosition>(MaxDuplicates, NullPosition);
  const inline static CellStoragePosition NullCellStoragePosition =
      CellStoragePosition(PreNullCellStoragePosition.begin(), PreNullCellStoragePosition.end());
  // TODO(David): initialize with NullPosition
  std::vector<CellStoragePosition> data = std::vector<CellStoragePosition>(NullCellStoragePosition);

  public:
  StorageBackmap(std::size_t size, LtsTree& storage)
      : storage(std::ref(storage)), data(NullCellStoragePosition, size) {}

  [[nodiscard]] std::size_t storageLookup(std::size_t id, std::size_t duplicate = 0) const {
    return storagePositionLookup(id, duplicate).second;
  }

  [[nodiscard]] StoragePosition storagePositionLookup(std::size_t id,
                                                      std::size_t duplicate = 0) const {
    assert(duplicate < MaxDuplicates);
    return data[id][duplicate];
  }

  template <typename TRef>
  std::size_t addElement(std::size_t color,
                         const TRef* zeroPosition,
                         const TRef* layerData,
                         int cell,
                         std::size_t index) {
    auto position = static_cast<std::size_t>((layerData + index) - zeroPosition) / sizeof(TRef);
    auto storagePosition = StoragePosition(color, position);
    std::size_t j = 0;
    for (j = 0; j < MaxDuplicates; ++j) {
      if (data[cell][j] == NullPosition) {
        data[cell][j] = storagePosition;
        break;
      }
    }
    assert(j < MaxDuplicates);
    return j;
  }

  template <typename TRef>
  void addLayer(std::size_t color,
                const TRef* zeroPosition,
                const TRef* layerData,
                const std::vector<std::size_t>& cells) {
    for (std::size_t i = 0; i < cells.size(); ++i) {
      addElement<TRef>(color, zeroPosition, layerData, cells[i], i);
    }
  }

  private:
  std::reference_wrapper<LTSTree> storage;
};

using ClusterBackmap = StorageBackmap<4, LTSTree>;

using DynrupBackmap = StorageBackmap<1, LTSTree>;

} // namespace seissol::initializer
#endif // SEISSOL_SRC_MEMORY_TREE_BACKMAP_H_
