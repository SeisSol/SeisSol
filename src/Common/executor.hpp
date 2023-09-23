#pragma once

namespace seissol {

enum class Executor { Host, Device };

constexpr bool executorEnabled(Executor executor) {
  // right now, SeisSol is only compiled either for CPU or for GPU---but never for both
#ifdef ACL_DEVICE
  return executor == Executor::Device;
#else
  return executor == Executor::Host;
#endif
}

} // namespace seissol
