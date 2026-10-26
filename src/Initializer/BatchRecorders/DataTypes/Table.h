// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_TABLE_H_
#define SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_TABLE_H_

#include "Condition.h"
#include "EncodedConstants.h"
#include "Memory/MemoryAllocator.h"

#include <array>
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

  std::vector<Type> getHostData() { return hostVector_; }
  [[nodiscard]] const std::vector<Type>& getHostData() const { return hostVector_; }

  size_t getSize() { return hostVector_.size(); }

  private:
  std::vector<Type> hostVector_;
  memory::MemkindArray<Type> deviceArray_;
};

/// The variables `KeyType::Id` of a batch. A layer records its batches in its configuration, so a
/// variable holds what that configuration holds, e.g. pointers to its reals: `set` takes the type
/// of the entries from the vector, and `get` names it again.
template <typename KeyType>
struct GenericTable {
  using VariableIdType = typename KeyType::Id;

  public:
  GenericTable() = default;

  template <typename DataType>
  void set(VariableIdType id, std::vector<DataType>& data) {
    content_[*id] = std::make_shared<GenericTableEntry<DataType>>(data);
    types_[*id] = &typeid(DataType);
  }

  template <typename DataType>
  GenericTableEntry<DataType>* get(VariableIdType id) {
    assert((types_.at(*id) == nullptr || *types_.at(*id) == typeid(DataType)) &&
           "a batch variable is read as another type than it was recorded as");
    return static_cast<GenericTableEntry<DataType>*>(content_.at(*id).get());
  }

  template <typename DataType>
  [[nodiscard]] const GenericTableEntry<DataType>* get(VariableIdType id) const {
    assert((types_.at(*id) == nullptr || *types_.at(*id) == typeid(DataType)) &&
           "a batch variable is read as another type than it was recorded as");
    return static_cast<const GenericTableEntry<DataType>*>(content_.at(*id).get());
  }

  private:
  std::array<std::shared_ptr<void>, *VariableIdType::Count> content_{};
  std::array<const std::type_info*, *VariableIdType::Count> types_{};
};

using PointersToRealsTable = GenericTable<inner_keys::Wp>;
using DrPointersToRealsTable = GenericTable<inner_keys::Dr>;
using IndicesTable = GenericTable<inner_keys::Indices>;

} // namespace seissol::recording

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_TABLE_H_
