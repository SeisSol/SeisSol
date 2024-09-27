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
 * Tree for managing lts data.
 **/

#ifndef INITIALIZER_TREE_LTSINTERNALNODE_HPP_
#define INITIALIZER_TREE_LTSINTERNALNODE_HPP_

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

    inline void nextLeaf() {
      do {
        ConstIterator::operator++();
      } while (*this != m_end && !m_node->isLeaf());
    }

    // m_node must point to a leaf or NULL
    inline void skipMaskedLayer() {
      while (*this != m_end && operator*().isMasked(m_layerMask)) {
        nextLeaf();
      }
    }

public:
    ConstLeafIterator(const ConstIterator& end) : ConstIterator(end) {}
    ConstLeafIterator(const ConstIterator& begin, const ConstIterator& end, LayerMask layerMask)
        : ConstIterator(begin), m_end(end), m_layerMask(layerMask) {}

    inline ConstLeafIterator& operator++() {
      nextLeaf();
      skipMaskedLayer();
      return *this;
    }

    inline const Layer& operator*() const { return *static_cast<const Layer*>(m_node); }

    inline const Layer* operator->() const { return static_cast<const Layer*>(m_node); }
  };

  inline LeafIterator beginLeaf(LayerMask layerMask = LayerMask()) {
    LeafIterator it = LeafIterator(begin(), end(), layerMask);
    it.skipMaskedLayer();
    return it;
  }

  inline LeafIterator endLeaf() { return LeafIterator(end()); }

  [[nodiscard]] inline ConstLeafIterator beginLeaf(LayerMask layerMask = LayerMask()) const {
    ConstLeafIterator it = ConstLeafIterator(begin(), end(), layerMask);
    it.skipMaskedLayer();
    return it;
  }

  [[nodiscard]] inline ConstLeafIterator endLeaf() const { return ConstLeafIterator(end()); }

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

    inline LeafIterator begin() { return node.beginLeaf(mask); }

    inline LeafIterator end() { return node.endLeaf(); }
  };

  class LeafIteratorWrapperConst {
private:
    const LTSInternalNode& node;
    LayerMask mask;

public:
    LeafIteratorWrapperConst(const LTSInternalNode& node, LayerMask mask)
        : node(node), mask(mask) {}

    [[nodiscard]] inline ConstLeafIterator begin() const { return node.beginLeaf(mask); }

    [[nodiscard]] inline ConstLeafIterator end() const { return node.endLeaf(); }
  };

  inline LeafIteratorWrapper leaves(LayerMask mask = LayerMask()) {
    return LeafIteratorWrapper(*this, mask);
  }

  [[nodiscard]] inline LeafIteratorWrapperConst leaves(LayerMask mask = LayerMask()) const {
    return LeafIteratorWrapperConst(*this, mask);
  }

  template <typename F>
  inline void iterateLeaves(F&& leafFunction, LayerMask mask = LayerMask()) {
    for (auto& leaf : leaves(mask)) {
      leafFunction(leaf);
    }
  }

  template <typename F>
  inline void iterateLeaves(F&& leafFunction, LayerMask mask = LayerMask()) const {
    for (const auto& leaf : leaves(mask)) {
      leafFunction(leaf);
    }
  }
};

} // namespace seissol::initializer

#endif
