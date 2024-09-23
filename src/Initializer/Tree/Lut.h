// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 */

#ifndef SEISSOL_SRC_INITIALIZER_TREE_LUT_H_
#define SEISSOL_SRC_INITIALIZER_TREE_LUT_H_

#include "LTSTree.h"
#include "Layer.h"

namespace seissol {
namespace initializer {
class Lut;
} // namespace initializer
} // namespace seissol

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
    unsigned* ltsToMesh;
    /** meshToLts[0][meshId] always returns a valid ltsId.
     * meshToLts[1..3][meshId] might be invalid (== std::numeric_limits<unsigned>::max())
     * and contains the ltsIds of duplicated cells.
     */
    unsigned* meshToLts[MaxDuplicates];
    /** Contains meshIds where any of meshToLts[1..3][meshId] is valid. */
    unsigned* duplicatedMeshIds;
    /** Size of duplicatedMeshIds. */
    unsigned numberOfDuplicatedMeshIds;

    LutsForMask();
    ~LutsForMask();

    void createLut(LayerMask mask,
                   LTSTree* ltsTree,
                   const unsigned* globalLtsToMesh,
                   unsigned numberOfMeshIds);
  };

  LutsForMask maskedLuts[1 << NumLayers];
  LTSTree* m_ltsTree;
  unsigned* m_meshToClusters;
  std::vector<LayerType> m_meshToLayer;

  public:
  Lut();
  ~Lut();

  void createLuts(LTSTree* ltsTree, unsigned* ltsToMesh, unsigned numberOfMeshIds);

  inline unsigned meshId(LayerMask mask, unsigned ltsId) const {
    return maskedLuts[mask.to_ulong()].ltsToMesh[ltsId];
  }

  inline unsigned* getLtsToMeshLut(LayerMask mask) const {
    return maskedLuts[mask.to_ulong()].ltsToMesh;
  }

  inline unsigned ltsId(LayerMask mask, unsigned meshId, unsigned duplicate = 0) const {
    assert(duplicate < MaxDuplicates);
    return maskedLuts[mask.to_ulong()].meshToLts[duplicate][meshId];
  }

  inline unsigned* const (&getMeshToLtsLut(LayerMask mask) const)[MaxDuplicates] {
    return maskedLuts[mask.to_ulong()].meshToLts;
  }

  inline unsigned* getDuplicatedMeshIds(LayerMask mask) const {
    return maskedLuts[mask.to_ulong()].duplicatedMeshIds;
  }

  inline unsigned getNumberOfDuplicatedMeshIds(LayerMask mask) const {
    return maskedLuts[mask.to_ulong()].numberOfDuplicatedMeshIds;
  }

  inline unsigned cluster(unsigned meshId) const { return m_meshToClusters[meshId]; }

  inline LayerType layer(unsigned meshId) const { return m_meshToLayer[meshId]; }

  inline unsigned* getMeshToClusterLut() const { return m_meshToClusters; }

  template <typename T>
  T& lookup(const Variable<T>& handle,
            unsigned meshId,
            AllocationPlace place = AllocationPlace::Host) const {
    return m_ltsTree->var(handle, place)[ltsId(handle.mask, meshId) * handle.count];
  }
};

#endif // SEISSOL_SRC_INITIALIZER_TREE_LUT_H_
