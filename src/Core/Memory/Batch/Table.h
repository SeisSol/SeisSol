// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_MEMORY_BATCH_TABLE_H_
#define SEISSOL_SRC_CORE_MEMORY_BATCH_TABLE_H_

#include "Condition.h"
#include "Core/Memory/MemoryAllocator.h"
#include "EncodedConstants.h"

#include <cassert>
#include <memory>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <utility>
#include <vector>

namespace seissol::recording {

template <typename Type>
class GenericTableEntry {
  using PointerType = Type*;

  public:
  explicit GenericTableEntry(const std::vector<Type>& userVector)
      : hostVector_(userVector), deviceArray_(userVector, memory::Memkind::DeviceGlobalMemory) {}

  PointerType getDeviceDataPtr() {
    assert(deviceArray_.data() != nullptr && "requested batch has not been recorded");
    return deviceArray_.data();
  }

  // a reference also for non-const entries: e.g. a host function enqueued on a stream still reads
  // the data after its caller has returned
  [[nodiscard]] const std::vector<Type>& getHostData() const { return hostVector_; }

  size_t getSize() { return hostVector_.size(); }

  private:
  std::vector<Type> hostVector_;
  memory::MemkindArray<Type> deviceArray_;
};

/// The variables of a batch (see BatchVariable). A layer records its batches in its configuration,
/// so a variable holds what that configuration holds, e.g. pointers to its reals: `set` takes the
/// type of the entries from the vector, and `get` names it again.
class BatchTable {
  public:
  template <typename Tag, typename DataType>
  void set(BatchVariable<Tag> /*variable*/, std::vector<DataType>& data) {
    const auto id = BatchVariable<Tag>::id();
    if (id >= content_.size()) {
      content_.resize(id + 1);
      types_.resize(id + 1, nullptr);
    }
    content_[id] = std::make_shared<GenericTableEntry<DataType>>(data);
    types_[id] = &typeid(DataType);
  }

  /// The entries of `variable`; null if the batch has none.
  template <typename DataType, typename Tag>
  GenericTableEntry<DataType>* get(BatchVariable<Tag> /*variable*/) {
    const auto id = BatchVariable<Tag>::id();
    if (id >= content_.size()) {
      return nullptr;
    }
    assert((types_[id] == nullptr || *types_[id] == typeid(DataType)) &&
           "a batch variable is read as another type than it was recorded as");
    return static_cast<GenericTableEntry<DataType>*>(content_[id].get());
  }

  /// The entries of `variable`; null if the batch has none.
  template <typename DataType, typename Tag>
  [[nodiscard]] const GenericTableEntry<DataType>* get(BatchVariable<Tag> /*variable*/) const {
    const auto id = BatchVariable<Tag>::id();
    if (id >= content_.size()) {
      return nullptr;
    }
    assert((types_[id] == nullptr || *types_[id] == typeid(DataType)) &&
           "a batch variable is read as another type than it was recorded as");
    return static_cast<const GenericTableEntry<DataType>*>(content_[id].get());
  }

  private:
  // by the id of the variable
  std::vector<std::shared_ptr<void>> content_;
  std::vector<const std::type_info*> types_;
};

} // namespace seissol::recording

#endif // SEISSOL_SRC_CORE_MEMORY_BATCH_TABLE_H_
