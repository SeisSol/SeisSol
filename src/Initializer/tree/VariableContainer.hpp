#ifndef INITIALIZER_TREE_VARIABLE_CONTAINER_HPP_
#define INITIALIZER_TREE_VARIABLE_CONTAINER_HPP_

namespace seissol::initializers {
class LTSTree;

struct LTSVariableContainer {
  virtual void addTo(LTSTree& tree) = 0;
};
} // namespace seissol::initializers

#endif
