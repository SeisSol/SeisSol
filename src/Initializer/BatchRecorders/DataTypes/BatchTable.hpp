#ifndef SEISSOL_POINTERSTABLE_HPP
#define SEISSOL_POINTERSTABLE_HPP

#ifdef ACL_DEVICE

#include "Condition.hpp"
#include "EncodedConstants.hpp"
#include <device.h>
#include <string>
#include <unordered_map>
#include <utility>
#include <array>
#include <vector>

namespace seissol::initializers::recording {

class BatchPointers {
public:
  explicit BatchPointers(std::vector<real *> collectedPointers)
      : pointers(std::move(collectedPointers)), devicePtrs(nullptr) {
    if (!pointers.empty()) {
      devicePtrs = (real **)device.api->allocGlobMem(pointers.size() * sizeof(real *));
      device.api->copyTo(devicePtrs, pointers.data(), pointers.size() * sizeof(real *));
    }
  }

  BatchPointers(const BatchPointers &other) : pointers(other.pointers), devicePtrs(nullptr) {
    if (!pointers.empty()) {
      if (other.devicePtrs != nullptr) {
        devicePtrs = (real **)device.api->allocGlobMem(other.pointers.size() * sizeof(real *));
        device.api->copyBetween(devicePtrs, other.devicePtrs, other.pointers.size() * sizeof(real *));
      }
    }
  }

  BatchPointers &operator=(const BatchPointers &Other) = delete;

  virtual ~BatchPointers() {
    if (devicePtrs != nullptr) {
      device.api->freeMem(devicePtrs);
      devicePtrs = nullptr;
    }
  }

  real **getPointers() {
    assert(devicePtrs != nullptr && "requested batch has not been recorded");
    return devicePtrs;
  }
  size_t getSize() {
    return pointers.size();
  }

private:
  std::vector<real *> pointers{};
  real **devicePtrs{nullptr};
  device::DeviceInstance &device = device::DeviceInstance::getInstance();
};

/**
 * This class may seem redundant. But it provides strong guarantee of
 * zero initialization of std::array. Note, there are some circumstances
 * when it is not zero-initialized
 * */
struct BatchTable {
public:
  BatchTable() {
    for (auto &item : content) {
      item = nullptr;
    }
  }
  ~BatchTable() {
    for (auto &item : content) {
      delete item;
    }
  }

  std::array<BatchPointers *, *EntityId::Count> content{};
};

} // namespace seissol::initializers::recording

#else  // ACL_DEVICE
namespace seissol::initializers::recording {
// Provide a dummy implementation for a pure CPU execution
struct BatchTable {};
} // namespace seissol::initializers::recording
#endif // ACL_DEVICE

#endif // SEISSOL_POINTERSTABLE_HPP
