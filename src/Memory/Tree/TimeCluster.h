// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_TREE_TIMECLUSTER_H_
#define SEISSOL_SRC_INITIALIZER_TREE_TIMECLUSTER_H_

#include "LTSInternalNode.h"
#include "Log2.h"

namespace seissol::initializer {

class TimeCluster : public LTSInternalNode {
  public:
  TimeCluster() {
    setChildren<Layer>(3);
    child<Ghost>().setLayerType(Ghost);
    child<Copy>().setLayerType(Copy);
    child<Interior>().setLayerType(Interior);
  }

  template <enum LayerType LAYER>
  Layer& child() {
    return *static_cast<Layer*>(m_children[Log2<LAYER>::Result].get());
  }

  Layer& child(LayerType type) {
    switch (type) {
    case Ghost:
      return child<Ghost>();
    case Copy:
      return child<Copy>();
    case Interior:
      return child<Interior>();
    default:
      throw;
    }
  }
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_TREE_TIMECLUSTER_H_
