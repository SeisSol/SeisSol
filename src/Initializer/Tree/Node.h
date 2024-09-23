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

#ifndef SEISSOL_SRC_INITIALIZER_TREE_NODE_H_
#define SEISSOL_SRC_INITIALIZER_TREE_NODE_H_

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

  class Iterator {
public:
    // NOLINTNEXTLINE
    using iterator_category = std::input_iterator_tag;
    // NOLINTNEXTLINE
    using value_type = Node;
    // NOLINTNEXTLINE
    using difference_type = size_t;
    // NOLINTNEXTLINE
    using pointer = Node*;
    // NOLINTNEXTLINE
    using reference = Node&;

    Iterator() : m_node(NULL) {}
    Iterator(Node* node) : m_node(node) {}

    inline Iterator& operator++() {
      m_node = m_node->m_next;
      return *this;
    }

    inline reference operator*() { return *m_node; }

    inline pointer operator->() { return m_node; }

    inline bool operator==(const Iterator& other) const { return other.m_node == m_node; }

    inline bool operator!=(const Iterator& other) const { return other.m_node != m_node; }

protected:
    value_type* m_node;
  };

  inline Iterator begin() {
    Node* start = this;
    while (!start->isLeaf()) {
      start = start->m_children[0].get();
    }
    assert(start == this || start->m_next != NULL);
    return Iterator(start);
  }

  inline Iterator end() {
    // The current node is the last one in a post-order traversal.
    // Hence, end() points to the node after the last one.
    return Iterator(this->m_next);
  }
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_TREE_NODE_H_
