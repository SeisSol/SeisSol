#pragma once

#include <Common/Executor.h>
#include <Initializer/Typedefs.h>
#include <Parallel/Runtime/Stream.h>

namespace seissol::kernels {

class Kernel {
  public:
  virtual void setGlobalData(const CompoundGlobalData& global) {}
  virtual ~Kernel() = default;

#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif
};

} // namespace seissol::kernels
