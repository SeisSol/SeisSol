#pragma once

namespace seissol {

enum class Executor { Host, Device };

constexpr bool executorEnabled(Executor executor) {
#ifdef ACL_DEVICE
  return true;
#else
  return executor == Executor::Host;
#endif
}

} // namespace seissol
