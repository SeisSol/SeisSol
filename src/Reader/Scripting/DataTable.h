// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_SCRIPTING_DATATABLE_H_
#define SEISSOL_SRC_READER_SCRIPTING_DATATABLE_H_

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>
namespace seissol::reader::scripting {

enum class DataType : std::uint8_t { F32, F64, I32, I64 };

template <typename T>
struct DataTypeTraits {
  static_assert(sizeof(T) == 0, "Unsupported type for scripting.");
};

template <>
struct DataTypeTraits<float> {
  static constexpr DataType Type = DataType::F32;
};
template <>
struct DataTypeTraits<double> {
  static constexpr DataType Type = DataType::F64;
};
template <>
struct DataTypeTraits<std::int32_t> {
  static constexpr DataType Type = DataType::I32;
};
template <>
struct DataTypeTraits<std::int64_t> {
  static constexpr DataType Type = DataType::I64;
};

enum class Direction : std::uint8_t { In, Out, InOut };

struct DataEntry {
  std::string name;
  Direction direction;
  DataType datatype{DataType::F64};
  std::function<void(std::size_t, void*)> accessor;
  std::function<void(std::size_t, const void*)> setter;

  template <typename T>
  T getValue(std::size_t index) const {
    assert(DataTypeTraits<T>::Type == datatype);

    T out{};
    accessor(index, &out);
    return out;
  }

  template <typename T>
  void getValueRange(T* storage, std::size_t start, std::size_t size) const {
    for (std::size_t i = 0; i < size; ++i) {
      storage[i] = getValue<T>(start + i);
    }
  }
};

class DataTable {
  public:
  explicit DataTable(std::size_t numPoints) : numPoints_(numPoints) {}

  // View-on-existing-storage
  template <typename T>
  void bindView(
      std::string name, Direction dir, T* base, std::size_t stride = 1, std::size_t offset = 0) {
    const auto accessor = [=](std::size_t idx, void* out) {
      auto* outC = reinterpret_cast<T*>(out);
      *outC = base[idx * stride + offset];
    };
    const auto setter = [=](std::size_t idx, const void* in) {
      auto* inC = reinterpret_cast<const T*>(in);
      base[idx * stride + offset] = *inC;
    };

    dataEntries_.emplace_back(
        DataEntry{std::move(name), dir, DataTypeTraits<T>::Type, accessor, setter});
  }

  // View-on-existing-struct
  template <typename S, typename T>
  void bindMemberView(std::string name, Direction dir, S* base, typename S::T* member) {
    const auto accessor = [=](std::size_t idx, void* out) {
      auto* outC = reinterpret_cast<T*>(out);
      *outC = base[idx].*member;
    };
    const auto setter = [=](std::size_t idx, const void* in) {
      auto* inC = reinterpret_cast<const T*>(in);
      base[idx].*member = *inC;
    };

    dataEntries_.emplace_back(
        DataEntry{std::move(name), dir, DataTypeTraits<T>::Type, accessor, setter});
  }

  // Lazy/computed (only called when reading)
  // signature (let's enforce it only once C++20 drops): (std::size_t index) -> returnType (i.e.
  // float/int)
  template <typename F>
  void bindComputed(std::string name, F&& fn) {
    using ReturnT = std::invoke_result_t<F, std::size_t>;

    const auto accessor = [=, ffn = std::forward<F>(fn)](std::size_t idx, void* out) {
      auto* outC = reinterpret_cast<ReturnT*>(out);
      *outC = ffn(idx);
    };

    dataEntries_.emplace_back(DataEntry{
        std::move(name), Direction::In, DataTypeTraits<ReturnT>::Type, accessor, nullptr});
  }

  [[nodiscard]] std::size_t numPoints() const { return numPoints_; }

  [[nodiscard]] const std::vector<DataEntry>& dataEntries() const { return dataEntries_; }

  private:
  std::size_t numPoints_;
  std::vector<DataEntry> dataEntries_;
};

} // namespace seissol::reader::scripting
#endif // SEISSOL_SRC_READER_SCRIPTING_DATATABLE_H_
