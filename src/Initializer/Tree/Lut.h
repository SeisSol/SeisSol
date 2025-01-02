/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Handles mapping between mesh and cells.
 **/

#ifndef INITIALIZER_TREE_LUT_HPP_
#define INITIALIZER_TREE_LUT_HPP_

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

#endif
