// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_TABLE_H_
#define SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_TABLE_H_

#ifdef ACL_DEVICE

#include "Condition.h"
#include "EncodedConstants.h"

#include <Device/device.h>
#include <array>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace seissol::recording {

template <typename Type>
class GenericTableEntry {
  using PointerType = Type*;

  public:
  explicit GenericTableEntry(std::vector<Type> userVector)
      : hostVector_(std::move(userVector)), deviceDataPtr_(nullptr) {
    if (!hostVector_.empty()) {
      deviceDataPtr_ =
          static_cast<PointerType>(device_.api->allocGlobMem(hostVector_.size() * sizeof(Type)));
      device_.api->copyTo(deviceDataPtr_, hostVector_.data(), hostVector_.size() * sizeof(Type));
    }
  }

  GenericTableEntry(const GenericTableEntry& other)
      : hostVector_(other.hostVector_), deviceDataPtr_(nullptr) {
    if (!hostVector_.empty()) {
      if (other.deviceDataPtr_ != nullptr) {
        deviceDataPtr_ = static_cast<PointerType>(
            device_.api->allocGlobMem(other.hostVector_.size() * sizeof(Type)));
        device_.api->copyBetween(
            deviceDataPtr_, other.deviceDataPtr_, other.hostVector_.size() * sizeof(Type));
      }
    }
  }

  GenericTableEntry& operator=(const GenericTableEntry& other) = delete;

  virtual ~GenericTableEntry() {
    if (deviceDataPtr_ != nullptr) {
      device_.api->freeGlobMem(deviceDataPtr_);
      deviceDataPtr_ = nullptr;
    }
  }

  PointerType getDeviceDataPtr() {
    assert(deviceDataPtr != nullptr && "requested batch has not been recorded");
    return deviceDataPtr_;
  }

  std::vector<Type> getHostData() { return hostVector_; }
  [[nodiscard]] const std::vector<Type>& getHostData() const { return hostVector_; }

  size_t getSize() { return hostVector_.size(); }

  private:
  std::vector<Type> hostVector_{};
  PointerType deviceDataPtr_{nullptr};
  device::DeviceInstance& device_ = device::DeviceInstance::getInstance();
};

template <typename KeyType>
struct GenericTable {
  using VariableIdType = typename KeyType::Id;
  using DataType = typename KeyType::DataType;

  public:
  GenericTable() = default;

  void set(VariableIdType id, std::vector<DataType>& data) {
    content_[*id] = std::make_shared<GenericTableEntry<DataType>>(data);
  }

  auto get(VariableIdType id) { return content_.at(*id).get(); }

  [[nodiscard]] auto get(VariableIdType id) const { return content_.at(*id).get(); }

  private:
  std::array<std::shared_ptr<GenericTableEntry<DataType>>, *VariableIdType::Count> content_{};
};

using PointersToRealsTable = GenericTable<inner_keys::Wp>;
using DrPointersToRealsTable = GenericTable<inner_keys::Dr>;
using MaterialTable = GenericTable<inner_keys::Material>;
using IndicesTable = GenericTable<inner_keys::Indices>;

} // namespace seissol::recording

#else  // ACL_DEVICE
namespace seissol::recording {
// Provide a dummy implementations for a pure CPU execution
struct PointersToRealsTable {};
struct DrPointersToRealsTable {};
struct MaterialTable {};
struct IndicesTable {};
} // namespace seissol::recording
#endif // ACL_DEVICE

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_TABLE_H_
