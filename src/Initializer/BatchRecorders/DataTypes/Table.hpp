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
template <typename KeyType>
struct GenericTable {
  using VariableIdType = typename KeyType::Id;
  using DataType = typename KeyType::DataType;

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

  void set(VariableIdType id, std::vector<DataType>& data) {
    content[*id] = new GenericTableEntry<DataType>(data);
  }

  auto get(VariableIdType id) { return content[*id]; }

  private:
  std::array<GenericTableEntry<DataType>*, *VariableIdType::Count> content{};
};

using PointersToRealsTable = GenericTable<inner_keys::Wp>;
using DrPointersToRealsTable = GenericTable<inner_keys::Dr>;
using MaterialTable = GenericTable<inner_keys::Material>;
using IndicesTable = GenericTable<inner_keys::Indices>;

} // namespace seissol::initializers::recording

#else  // ACL_DEVICE
namespace seissol::initializers::recording {
// Provide a dummy implementations for a pure CPU execution
struct PointersToRealsTable {};
struct DrPointersToRealsTable {};
struct MaterialTable {};
struct IndicesTable {};
} // namespace seissol::initializers::recording
#endif // ACL_DEVICE

#endif // SEISSOL_POINTERSTABLE_HPP
