#pragma once

#include <cstddef>
namespace seissol {

class OpenMP {
  public:
  static bool enabled();
  static std::size_t threadId();
  static std::size_t threadCount();
};

} // namespace seissol
