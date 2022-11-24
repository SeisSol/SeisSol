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

template <typename Type>
class GenericTableEntry {
  using PointerType = Type*;

  public:
  explicit GenericTableEntry(std::vector<Type> userVector)
      : hostVector(std::move(userVector)), deviceDataPtr(nullptr) {
    if (!hostVector.empty()) {
      deviceDataPtr =
          static_cast<PointerType>(device.api->allocGlobMem(hostVector.size() * sizeof(Type)));
      device.api->copyTo(deviceDataPtr, hostVector.data(), hostVector.size() * sizeof(Type));
    }
  }

  GenericTableEntry(const GenericTableEntry& other)
      : hostVector(other.hostVector), deviceDataPtr(nullptr) {
    if (!hostVector.empty()) {
      if (other.devicePtrs != nullptr) {
        deviceDataPtr = static_cast<PointerType>(
            device.api->allocGlobMem(other.pointers.size() * sizeof(Type)));
        device.api->copyBetween(
            deviceDataPtr, other.devicePtrs, other.pointers.size() * sizeof(Type));
      }
    }
  }

  GenericTableEntry& operator=(const GenericTableEntry& other) = delete;

  virtual ~GenericTableEntry() {
    if (deviceDataPtr != nullptr) {
      device.api->freeMem(deviceDataPtr);
      deviceDataPtr = nullptr;
    }
  }

  PointerType getDeviceDataPtr() {
    assert(devicePtrs != nullptr && "requested batch has not been recorded");
    return deviceDataPtr;
  }

  std::vector<Type> getHostData() { return hostVector; }

  const std::vector<Type>& getHostData() const { return hostVector; }

  size_t getSize() { return hostVector.size(); }

  private:
  std::vector<Type> hostVector{};
  PointerType deviceDataPtr{nullptr};
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
};

/**
 * This class may seem redundant. But it provides strong guarantee of
 * zero initialization of std::array. Note, there are some circumstances
 * when it is not zero-initialized
 * */
template <typename T>
struct GenericTable {
  public:
  GenericTable() {
    for (auto& entry : content) {
      entry = nullptr;
    }
  }
  ~GenericTable() {
    for (auto& entry : content) {
      delete entry;
    }
  }

  void set(enum EntityId id, std::vector<T>& data) {
    content[*id] = new GenericTableEntry<T>(data);
  }

  auto get(enum EntityId id) { return content[*id]; }

  private:
  std::array<GenericTableEntry<T>*, *EntityId::Count> content{};
};

using PointersToRealsTable = GenericTable<real*>;
using IndicesTable = GenericTable<unsigned>;

} // namespace seissol::initializers::recording

#else  // ACL_DEVICE
namespace seissol::initializers::recording {
// Provide a dummy implementation for a pure CPU execution
struct PointersToRealsTable {};
struct IndicesTable {};
} // namespace seissol::initializers::recording
#endif // ACL_DEVICE

#endif // SEISSOL_POINTERSTABLE_HPP
