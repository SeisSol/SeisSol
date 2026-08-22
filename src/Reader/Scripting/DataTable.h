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
#include <memory>
#include <optional>
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

/// The address arithmetic behind a bound column, for consumers that cannot go
/// through the accessor.
///
/// ADDED (reported, Package 5). Two of them need it. A compiled kernel called
/// once per element must not pay a std::function per point -- measured at 0.4 us
/// per face just to BUILD the table, against 0.3 us to evaluate it, so the
/// indirection already costs more than the arithmetic. And a device kernel
/// cannot call a std::function at all, at any price.
///
/// Only the view forms can fill this in: bindView and bindMemberView are both
/// `base + index * byteStride + byteOffset` once the element type is erased,
/// and bindConstant is the same with a stride of zero. bindComputed cannot,
/// which is exactly what makes it the form that disqualifies a program from the
/// device path -- see makeKernel.
struct StridedView {
  void* base{nullptr}; // const-ness is carried by DataEntry::direction
  std::size_t byteStride{0};
  std::size_t byteOffset{0};
  bool writable{false};
};

struct DataEntry {
  std::string name;
  Direction direction;
  DataType datatype{DataType::F64};
  std::function<void(std::size_t, void*)> accessor;
  std::function<void(std::size_t, const void*)> setter;
  /// Empty for bindComputed; set by every other bind form.
  std::optional<StridedView> view;

  template <typename T>
  [[nodiscard]] T getValue(std::size_t index) const {
    assert(DataTypeTraits<T>::Type == datatype);

    T out{};
    accessor(index, &out);
    return out;
  }

  template <typename T>
  void setValue(std::size_t index, T value) const {
    assert(DataTypeTraits<T>::Type == datatype);

    setter(index, &value);
  }

  template <typename T>
  [[nodiscard]] T getValueAs(std::size_t index) const {
    switch (datatype) {
    case DataType::F32:
      return getValue<float>(index);
    case DataType::F64:
      return getValue<double>(index);
    case DataType::I32:
      return getValue<int32_t>(index);
    case DataType::I64:
      return getValue<int64_t>(index);
    }
    throw;
  }

  template <typename T>
  void setValueAs(std::size_t index, T value) const {
    switch (datatype) {
    case DataType::F32:
      setValue<float>(index, value);
      return;
    case DataType::F64:
      setValue<double>(index, value);
      return;
    case DataType::I32:
      setValue<int32_t>(index, value);
      return;
    case DataType::I64:
      setValue<int64_t>(index, value);
      return;
    }
    throw;
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

    dataEntries_.emplace_back(DataEntry{std::move(name),
                                        dir,
                                        DataTypeTraits<T>::Type,
                                        accessor,
                                        setter,
                                        makeView(base, stride, offset, true)});
  }

  /// A value that is the same at every point.
  ///
  /// ADDED (reported, Package 5). Expressible today only as bindComputed, which
  /// is the one form with no address arithmetic behind it -- so a table that
  /// merely names a constant (`group`, `sim` in EasiBoundary::query) cannot
  /// reach a device kernel. As a stride-0 view it can, and it is the same shape
  /// Package 2 chose for folding script parameters in as Const.
  template <typename T>
  void bindConstant(std::string name, const T& value) {
    auto held = std::make_shared<T>(value);
    const auto accessor = [held](std::size_t, void* out) { *reinterpret_cast<T*>(out) = *held; };
    StridedView view;
    view.base = const_cast<void*>(static_cast<const void*>(held.get()));
    view.byteStride = 0;
    view.byteOffset = 0;
    view.writable = false;
    dataEntries_.emplace_back(DataEntry{
        std::move(name), Direction::In, DataTypeTraits<T>::Type, accessor, nullptr, view});
    constants_.push_back(std::move(held));
  }

  // View-on-existing-storage
  template <typename T>
  void bindViewConst(std::string name,
                     Direction dir,
                     const T* base,
                     std::size_t stride = 1,
                     std::size_t offset = 0) {
    const auto accessor = [=](std::size_t idx, void* out) {
      auto* outC = reinterpret_cast<T*>(out);
      *outC = base[idx * stride + offset];
    };

    dataEntries_.emplace_back(DataEntry{std::move(name),
                                        dir,
                                        DataTypeTraits<T>::Type,
                                        accessor,
                                        nullptr,
                                        makeView(base, stride, offset, false)});
  }

  // View-on-existing-struct
  template <typename S, typename T>
  void bindMemberView(std::string name, Direction dir, S* base, T S::* member) {
    const auto accessor = [=](std::size_t idx, void* out) {
      auto* outC = reinterpret_cast<T*>(out);
      *outC = base[idx].*member;
    };
    const auto setter = [=](std::size_t idx, const void* in) {
      auto* inC = reinterpret_cast<const T*>(in);
      base[idx].*member = *inC;
    };

    dataEntries_.emplace_back(DataEntry{std::move(name),
                                        dir,
                                        DataTypeTraits<T>::Type,
                                        accessor,
                                        setter,
                                        makeMemberView(base, member, true)});
  }

  // View-on-existing-struct
  template <typename S, typename T>
  void bindMemberViewConst(std::string name, Direction dir, const S* base, T S::* member) {
    const auto accessor = [=](std::size_t idx, void* out) {
      auto* outC = reinterpret_cast<T*>(out);
      *outC = base[idx].*member;
    };

    dataEntries_.emplace_back(DataEntry{std::move(name),
                                        dir,
                                        DataTypeTraits<T>::Type,
                                        accessor,
                                        nullptr,
                                        makeMemberView(base, member, false)});
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

    // No StridedView on purpose: there is no address arithmetic behind a
    // functor, so this is the one form a device kernel cannot consume.
    dataEntries_.emplace_back(DataEntry{std::move(name),
                                        Direction::In,
                                        DataTypeTraits<ReturnT>::Type,
                                        accessor,
                                        nullptr,
                                        std::nullopt});
  }

  [[nodiscard]] std::size_t numPoints() const { return numPoints_; }

  [[nodiscard]] const std::vector<DataEntry>& dataEntries() const { return dataEntries_; }

  private:
  /// Element-typed stride/offset erased to bytes, which is what both consumers
  /// of StridedView actually need and the only form a device kernel can use.
  template <typename T>
  static StridedView makeView(T* base, std::size_t stride, std::size_t offset, bool writable) {
    StridedView view;
    view.base = const_cast<void*>(static_cast<const void*>(base));
    view.byteStride = stride * sizeof(T);
    view.byteOffset = offset * sizeof(T);
    view.writable = writable;
    return view;
  }

  /// A member view is a strided view with the struct as the stride: base[i].m
  /// is base + i * sizeof(S) + offsetof(S, m). Computed from a null object
  /// rather than offsetof because the member is a pointer-to-member, and this
  /// is well defined for the standard-layout structs the ABI is used with.
  template <typename S, typename T>
  static StridedView makeMemberView(S* base, T S::* member, bool writable) {
    const auto* null = static_cast<const S*>(nullptr);
    StridedView view;
    view.base = const_cast<void*>(static_cast<const void*>(base));
    view.byteStride = sizeof(S);
    view.byteOffset = static_cast<std::size_t>(reinterpret_cast<const char*>(&(null->*member)) -
                                               reinterpret_cast<const char*>(null));
    view.writable = writable;
    return view;
  }

  std::vector<std::shared_ptr<void>> constants_;
  std::size_t numPoints_;
  std::vector<DataEntry> dataEntries_;
};

} // namespace seissol::reader::scripting
#endif // SEISSOL_SRC_READER_SCRIPTING_DATATABLE_H_
