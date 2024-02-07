#pragma once

#include <functional>
#include <utility>

namespace seissol::parallel {
template <typename GraphT>
class ComputeGraph {
  private:
  bool captured = false;
  GraphT implementation;

  public:
  template <typename F>
  void run(typename GraphT::StreamT& stream, F&& handler) {
    if (implementation.canCapture()) {
      if (!captured) {
        implementation.beginCapture(stream);
        std::invoke(std::forward<F>(handler), stream);
        implementation.endCapture(stream);
        captured = true;
      }
      implementation.runCapture(stream);
    } else {
      std::invoke(std::forward<F>(handler), stream);
    }
  }
};
} // namespace seissol::parallel
