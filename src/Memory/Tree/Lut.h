// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_TREE_LUT_H_
#define SEISSOL_SRC_INITIALIZER_TREE_LUT_H_

#include "LTSTree.h"
#include "Layer.h"

namespace seissol::initializer {
class Lut;
} // namespace seissol::initializer

class seissol::initializer::Lut {
  public:
  static const unsigned MaxDuplicates = 4;

  private:
  /** Creates lookup tables (lut) for a given layer mask.
   *  ltsIds are consecutive and are to be understood with respect to the mask.
   *  I.e. if a variable is stored on the copy and the interior layer, then
   *  no ids are assigned to cells on the ghost layer.
   *
   *  meshIds might be invalid (== std::numeric_limits<unsigned>::max())
   * */
  struct LutsForMask {
    /** ltsToMesh[ltsId] returns a meshId given a ltsId. */
    std::vector<unsigned> ltsToMesh;
    /** meshToLts[0][meshId] always returns a valid ltsId.
     * meshToLts[1..3][meshId] might be invalid (== std::numeric_limits<unsigned>::max())
     * and contains the ltsIds of duplicated cells.
     */
    std::array<std::vector<unsigned>, MaxDuplicates> meshToLts;
    /** Contains meshIds where any of meshToLts[1..3][meshId] is valid. */
    std::vector<unsigned> duplicatedMeshIds;

    LutsForMask() = default;

    void createLut(LayerMask mask,
                   LTSTree* ltsTree,
                   const unsigned* globalLtsToMesh,
                   unsigned numberOfMeshIds);
  };

  LutsForMask maskedLuts[1 << NumLayers];
  LTSTree* m_ltsTree{nullptr};
  std::vector<unsigned> m_meshToClusters;
  std::vector<LayerType> m_meshToLayer;

  public:
  Lut();

  void createLuts(LTSTree* ltsTree, unsigned* ltsToMesh, unsigned numberOfMeshIds);

  [[nodiscard]] unsigned meshId(LayerMask mask, unsigned ltsId) const {
    return maskedLuts[mask.to_ulong()].ltsToMesh[ltsId];
  }

  [[nodiscard]] auto getLtsToMeshLut(LayerMask mask) const {
    return maskedLuts[mask.to_ulong()].ltsToMesh.data();
  }

  [[nodiscard]] unsigned ltsId(LayerMask mask, unsigned meshId, unsigned duplicate = 0) const {
    assert(duplicate < MaxDuplicates);
    return maskedLuts[mask.to_ulong()].meshToLts[duplicate][meshId];
  }

  [[nodiscard]] const auto& getMeshToLtsLut(LayerMask mask) const {
    return maskedLuts[mask.to_ulong()].meshToLts;
  }

  [[nodiscard]] auto getDuplicatedMeshIds(LayerMask mask) const {
    return maskedLuts[mask.to_ulong()].duplicatedMeshIds.data();
  }

  [[nodiscard]] unsigned getNumberOfDuplicatedMeshIds(LayerMask mask) const {
    return maskedLuts[mask.to_ulong()].duplicatedMeshIds.size();
  }

  [[nodiscard]] unsigned cluster(unsigned meshId) const { return m_meshToClusters[meshId]; }

  [[nodiscard]] LayerType layer(unsigned meshId) const { return m_meshToLayer[meshId]; }

  [[nodiscard]] const auto& getMeshToClusterLut() const { return m_meshToClusters; }

  template <typename T>
  T& lookup(const Variable<T>& handle,
            unsigned meshId,
            AllocationPlace place = AllocationPlace::Host) const {
    return m_ltsTree->var(handle, place)[ltsId(handle.mask, meshId) * handle.count];
  }
};

#endif // SEISSOL_SRC_INITIALIZER_TREE_LUT_H_
