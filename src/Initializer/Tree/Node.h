// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

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
  Node* m_next{nullptr};

  void setPostOrderPointers(Node* previous = nullptr) {
    for (auto& child : m_children) {
      child->setPostOrderPointers(previous);
      previous = child.get();
    }
    if (previous != nullptr) {
      previous->m_next = this;
    }
  }

  public:
  Node() = default;
  virtual ~Node() = default;

  template <typename T>
  void setChildren(unsigned numChildren) {
    m_children.resize(numChildren);
    for (auto& child : m_children) {
      child = std::make_shared<T>();
    }
  }

  [[nodiscard]] bool isLeaf() const { return m_children.empty(); }

  [[nodiscard]] unsigned numChildren() const { return m_children.size(); }

  class Iterator {
public:
    // NOLINTNEXTLINE
    using iterator_category = std::input_iterator_tag;
    // NOLINTNEXTLINE
    using value_type = Node;
    // NOLINTNEXTLINE
    using difference_type = ssize_t;
    // NOLINTNEXTLINE
    using pointer = Node*;
    // NOLINTNEXTLINE
    using reference = Node&;

    Iterator() : m_node(nullptr) {}
    Iterator(Node* node) : m_node(node) {}

    Iterator& operator++() {
      m_node = m_node->m_next;
      return *this;
    }

    reference operator*() { return *m_node; }

    pointer operator->() { return m_node; }

    bool operator==(const Iterator& other) const { return other.m_node == m_node; }

    bool operator!=(const Iterator& other) const { return other.m_node != m_node; }

protected:
    value_type* m_node;
  };

  class ConstIterator {
public:
    // NOLINTNEXTLINE
    using iterator_category = std::input_iterator_tag;
    // NOLINTNEXTLINE
    using value_type = const Node;
    // NOLINTNEXTLINE
    using difference_type = ssize_t;
    // NOLINTNEXTLINE
    using pointer = const Node*;
    // NOLINTNEXTLINE
    using reference = const Node&;

    ConstIterator() : m_node(nullptr) {}
    ConstIterator(const Node* node) : m_node(node) {}

    ConstIterator& operator++() {
      m_node = m_node->m_next;
      return *this;
    }

    reference operator*() { return *m_node; }

    pointer operator->() { return m_node; }

    bool operator==(const ConstIterator& other) const { return other.m_node == m_node; }

    bool operator!=(const ConstIterator& other) const { return other.m_node != m_node; }

protected:
    const value_type* m_node;
  };

  Iterator begin() {
    Node* start = this;
    while (!start->isLeaf()) {
      start = start->m_children[0].get();
    }
    assert(start == this || start->m_next != nullptr);
    return {start};
  }

  Iterator end() {
    // The current node is the last one in a post-order traversal.
    // Hence, end() points to the node after the last one.
    return {this->m_next};
  }

  [[nodiscard]] ConstIterator begin() const {
    const Node* start = this;
    while (!start->isLeaf()) {
      start = start->m_children[0].get();
    }
    assert(start == this || start->m_next != nullptr);
    return {start};
  }

  [[nodiscard]] ConstIterator end() const {
    // The current node is the last one in a post-order traversal.
    // Hence, end() points to the node after the last one.
    return {this->m_next};
  }
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_TREE_NODE_H_
