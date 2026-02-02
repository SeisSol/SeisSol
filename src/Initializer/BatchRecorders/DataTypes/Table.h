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

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_TABLE_H_
