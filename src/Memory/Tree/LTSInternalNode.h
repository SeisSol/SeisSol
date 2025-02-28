// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_TREE_LTSINTERNALNODE_H_
#define SEISSOL_SRC_INITIALIZER_TREE_LTSINTERNALNODE_H_

#include "Layer.h"
#include "Node.h"

namespace seissol::initializer {

class LTSInternalNode : public Node {
  public:
  class LeafIterator : public Iterator {
    friend class LTSInternalNode;

public:
    // NOLINTNEXTLINE
    using iterator_category = std::input_iterator_tag;
    // NOLINTNEXTLINE
    using value_type = Layer;
    // NOLINTNEXTLINE
    using difference_type = ssize_t;
    // NOLINTNEXTLINE
    using pointer = Layer*;
    // NOLINTNEXTLINE
    using reference = Layer&;

private:
    Iterator m_end;
    LayerMask m_layerMask;

    void nextLeaf() {
      do {
        Iterator::operator++();
      } while (*this != m_end && !m_node->isLeaf());
    }

    // m_node must point to a leaf or NULL
    void skipMaskedLayer() {
      while (*this != m_end && operator*().isMasked(m_layerMask)) {
        nextLeaf();
      }
    }

public:
    LeafIterator(const Iterator& end) : Iterator(end) {}
    LeafIterator(const Iterator& begin, const Iterator& end, LayerMask layerMask)
        : Iterator(begin), m_end(end), m_layerMask(layerMask) {}

    LeafIterator& operator++() {
      nextLeaf();
      skipMaskedLayer();
      return *this;
    }

    Layer& operator*() { return *dynamic_cast<Layer*>(m_node); }

    Layer* operator->() { return dynamic_cast<Layer*>(m_node); }
  };

  class ConstLeafIterator : public ConstIterator {
    friend class LTSInternalNode;

public:
    // NOLINTNEXTLINE
    using iterator_category = std::input_iterator_tag;
    // NOLINTNEXTLINE
    using value_type = const Layer;
    // NOLINTNEXTLINE
    using difference_type = ssize_t;
    // NOLINTNEXTLINE
    using pointer = const Layer*;
    // NOLINTNEXTLINE
    using reference = const Layer&;

private:
    ConstIterator m_end;
    LayerMask m_layerMask;

    void nextLeaf() {
      do {
        ConstIterator::operator++();
      } while (*this != m_end && !m_node->isLeaf());
    }

    // m_node must point to a leaf or NULL
    void skipMaskedLayer() {
      while (*this != m_end && operator*().isMasked(m_layerMask)) {
        nextLeaf();
      }
    }

public:
    ConstLeafIterator(const ConstIterator& end) : ConstIterator(end) {}
    ConstLeafIterator(const ConstIterator& begin, const ConstIterator& end, LayerMask layerMask)
        : ConstIterator(begin), m_end(end), m_layerMask(layerMask) {}

    ConstLeafIterator& operator++() {
      nextLeaf();
      skipMaskedLayer();
      return *this;
    }

    const Layer& operator*() const { return *dynamic_cast<const Layer*>(m_node); }

    const Layer* operator->() const { return dynamic_cast<const Layer*>(m_node); }
  };

  LeafIterator beginLeaf(LayerMask layerMask = LayerMask()) {
    LeafIterator it = LeafIterator(begin(), end(), layerMask);
    it.skipMaskedLayer();
    return it;
  }

  LeafIterator endLeaf() { return LeafIterator(end()); }

  [[nodiscard]] ConstLeafIterator beginLeaf(LayerMask layerMask = LayerMask()) const {
    ConstLeafIterator it = ConstLeafIterator(begin(), end(), layerMask);
    it.skipMaskedLayer();
    return it;
  }

  [[nodiscard]] ConstLeafIterator endLeaf() const { return ConstLeafIterator(end()); }

  [[nodiscard]] unsigned getNumberOfCells(LayerMask layerMask = LayerMask()) const {
    unsigned numCells = 0;
    for (const auto& leaf : leaves(layerMask)) {
      numCells += leaf.getNumberOfCells();
    }
    return numCells;
  }

  class LeafIteratorWrapper {
private:
    LTSInternalNode& node;
    LayerMask mask;

public:
    LeafIteratorWrapper(LTSInternalNode& node, LayerMask mask) : node(node), mask(mask) {}

    LeafIterator begin() { return node.beginLeaf(mask); }

    LeafIterator end() { return node.endLeaf(); }
  };

  class LeafIteratorWrapperConst {
private:
    const LTSInternalNode& node;
    LayerMask mask;

public:
    LeafIteratorWrapperConst(const LTSInternalNode& node, LayerMask mask)
        : node(node), mask(mask) {}

    [[nodiscard]] ConstLeafIterator begin() const { return node.beginLeaf(mask); }

    [[nodiscard]] ConstLeafIterator end() const { return node.endLeaf(); }
  };

  LeafIteratorWrapper leaves(LayerMask mask = LayerMask()) {
    return LeafIteratorWrapper(*this, mask);
  }

  [[nodiscard]] LeafIteratorWrapperConst leaves(LayerMask mask = LayerMask()) const {
    return LeafIteratorWrapperConst(*this, mask);
  }

  template <typename F>
  void iterateLeaves(F leafFunction, LayerMask mask = LayerMask()) {
    for (auto& leaf : leaves(mask)) {
      leafFunction(leaf);
    }
  }

  template <typename F>
  void iterateLeaves(F leafFunction, LayerMask mask = LayerMask()) const {
    for (const auto& leaf : leaves(mask)) {
      leafFunction(leaf);
    }
  }
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_TREE_LTSINTERNALNODE_H_
