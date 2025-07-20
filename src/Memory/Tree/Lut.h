// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MEMORY_TREE_LUT_H_
#define SEISSOL_SRC_MEMORY_TREE_LUT_H_

#include "Layer.h"
#include "Memory/Descriptor/LTS.h"

namespace seissol::initializer {
class Lut;
} // namespace seissol::initializer

class seissol::initializer::Lut {
  public:
  static const std::size_t MaxDuplicates = 4;

  private:
  /** Creates lookup tables (lut) for a given layer mask.
   *  ltsIds are consecutive and are to be understood with respect to the mask.
   *  I.e. if a variable is stored on the copy and the interior layer, then
   *  no ids are assigned to cells on the ghost layer.
   *
   *  meshIds might be invalid (== std::numeric_limits<std::size_t>::max())
   * */
  struct LutsForMask {
    /** ltsToMesh[ltsId] returns a meshId given a ltsId. */
    std::vector<std::size_t> ltsToMesh;
    /** meshToLts[0][meshId] always returns a valid ltsId.
     * meshToLts[1..3][meshId] might be invalid (== std::numeric_limits<std::size_t>::max())
     * and contains the ltsIds of duplicated cells.
     */
    std::array<std::vector<std::size_t>, MaxDuplicates> meshToLts;
    /** Contains meshIds where any of meshToLts[1..3][meshId] is valid. */
    std::vector<std::size_t> duplicatedMeshIds;

    LutsForMask() = default;

    void createLut(LayerMask mask,
                   LTS::Tree* ltsTree,
                   const std::size_t* globalLtsToMesh,
                   std::size_t numberOfMeshIds);
  };

  LutsForMask maskedLuts[1U << NumLayers];
  LTS::Tree* m_ltsTree{nullptr};
  std::vector<std::size_t> m_meshToClusters;
  std::vector<LayerType> m_meshToLayer;

  public:
  Lut();

  void createLuts(LTS::Tree* ltsTree, std::size_t* ltsToMesh, std::size_t numberOfMeshIds);

  [[nodiscard]] std::size_t meshId(LayerMask mask, std::size_t ltsId) const {
    return maskedLuts[mask.to_ulong()].ltsToMesh[ltsId];
  }

  [[nodiscard]] auto getLtsToMeshLut(LayerMask mask) const {
    return maskedLuts[mask.to_ulong()].ltsToMesh.data();
  }

  [[nodiscard]] std::size_t
      ltsId(LayerMask mask, std::size_t meshId, std::size_t duplicate = 0) const {
    assert(duplicate < MaxDuplicates);
    return maskedLuts[mask.to_ulong()].meshToLts[duplicate][meshId];
  }

  [[nodiscard]] const auto& getMeshToLtsLut(LayerMask mask) const {
    return maskedLuts[mask.to_ulong()].meshToLts;
  }

  [[nodiscard]] auto getDuplicatedMeshIds(LayerMask mask) const {
    return maskedLuts[mask.to_ulong()].duplicatedMeshIds.data();
  }

  [[nodiscard]] std::size_t getNumberOfDuplicatedMeshIds(LayerMask mask) const {
    return maskedLuts[mask.to_ulong()].duplicatedMeshIds.size();
  }

  [[nodiscard]] std::size_t cluster(std::size_t meshId) const { return m_meshToClusters[meshId]; }

  [[nodiscard]] LayerType layer(std::size_t meshId) const { return m_meshToLayer[meshId]; }

  [[nodiscard]] const auto& getMeshToClusterLut() const { return m_meshToClusters; }

  template <typename HandleT>
  typename HandleT::Type& lookup(const HandleT& handle,
                                 std::size_t meshId,
                                 AllocationPlace place = AllocationPlace::Host) const {
    return m_ltsTree->var(handle, place)[ltsId(m_ltsTree->info(handle).mask, meshId)];
  }

  template <typename StorageT>
  typename StorageT::Type& lookup(std::size_t meshId,
                                  AllocationPlace place = AllocationPlace::Host) const {
    return m_ltsTree->var<StorageT>(place)[ltsId(m_ltsTree->info<StorageT>().mask, meshId)];
  }
};

#endif // SEISSOL_SRC_MEMORY_TREE_LUT_H_
