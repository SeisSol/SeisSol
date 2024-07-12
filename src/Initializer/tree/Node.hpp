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
 * General purpose node for a data tree.
 **/

#ifndef INITIALIZER_TREE_NODE_HPP_
#define INITIALIZER_TREE_NODE_HPP_

#include <cassert>
#include <iterator>
#include <memory>
#include <vector>

namespace seissol::initializer {

class Node {
  protected:
  std::vector<std::shared_ptr<Node>> m_children;
  Node* m_next;

  void setPostOrderPointers(Node* previous = nullptr) {
    for (unsigned child = 0; child < m_children.size(); ++child) {
      m_children[child]->setPostOrderPointers(previous);
      previous = m_children[child].get();
    }
    if (previous != nullptr) {
      previous->m_next = this;
    }
  }

  public:
  Node() : m_children(), m_next(nullptr) {}
  virtual ~Node() = default;

  template <typename T>
  void setChildren(unsigned numChildren) {
    m_children.resize(numChildren);
    for (unsigned child = 0; child < m_children.size(); ++child) {
      m_children[child] = std::make_shared<T>();
    }
  }

  inline bool isLeaf() const { return m_children.empty(); }

  inline unsigned numChildren() const { return m_children.size(); }

  class iterator {
public:
    using iterator_category = std::input_iterator_tag;
    using value_type = Node;
    using difference_type = size_t;
    using pointer = Node*;
    using reference = Node&;

    iterator() : m_node(NULL) {}
    iterator(Node* node) : m_node(node) {}

    inline iterator& operator++() {
      m_node = m_node->m_next;
      return *this;
    }

    inline reference operator*() { return *m_node; }

    inline pointer operator->() { return m_node; }

    inline bool operator==(const iterator& other) const { return other.m_node == m_node; }

    inline bool operator!=(const iterator& other) const { return other.m_node != m_node; }

protected:
    value_type* m_node;
  };

  inline iterator begin() {
    Node* start = this;
    while (!start->isLeaf()) {
      start = start->m_children[0].get();
    }
    assert(start == this || start->m_next != NULL);
    return iterator(start);
  }

  inline iterator end() {
    // The current node is the last one in a post-order traversal.
    // Hence, end() points to the node after the last one.
    return iterator(this->m_next);
  }
};

} // namespace seissol::initializer

#endif
