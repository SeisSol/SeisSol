// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
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

#ifndef SEISSOL_SRC_INITIALIZER_TREE_LTSINTERNALNODE_H_
#define SEISSOL_SRC_INITIALIZER_TREE_LTSINTERNALNODE_H_

#include "Layer.h"
#include "Node.h"

namespace seissol {
namespace initializer {
class LTSInternalNode;
} // namespace initializer
} // namespace seissol

class seissol::initializer::LTSInternalNode : public seissol::initializer::Node {
  public:
  class LeafIterator : public Iterator {
    friend class LTSInternalNode;

private:
    Iterator m_end;
    LayerMask m_layerMask;

    inline void nextLeaf() {
      do {
        Iterator::operator++();
      } while (*this != m_end && !m_node->isLeaf());
    }

    // m_node must point to a leaf or NULL
    inline void skipMaskedLayer() {
      while (*this != m_end && operator*().isMasked(m_layerMask)) {
        nextLeaf();
      }
    }

public:
    LeafIterator(const Iterator& end) : Iterator(end) {}
    LeafIterator(const Iterator& begin, const Iterator& end, LayerMask layerMask)
        : Iterator(begin), m_end(end), m_layerMask(layerMask) {}

    inline LeafIterator& operator++() {
      nextLeaf();
      skipMaskedLayer();
      return *this;
    }

    inline Layer& operator*() { return *static_cast<Layer*>(m_node); }

    inline Layer* operator->() { return static_cast<Layer*>(m_node); }
  };

  inline LeafIterator beginLeaf(LayerMask layerMask = LayerMask()) {
    LeafIterator it = LeafIterator(begin(), end(), layerMask);
    it.skipMaskedLayer();
    return it;
  }

  inline LeafIterator endLeaf() { return LeafIterator(end()); }

  unsigned getNumberOfCells(LayerMask layerMask = LayerMask()) {
    unsigned numCells = 0;
    for (auto it = beginLeaf(layerMask); it != endLeaf(); ++it) {
      numCells += it->getNumberOfCells();
    }
    return numCells;
  }
};

#endif // SEISSOL_SRC_INITIALIZER_TREE_LTSINTERNALNODE_H_
