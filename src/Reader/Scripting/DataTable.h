// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_SCRIPTING_DATATABLE_H_
#define SEISSOL_SRC_READER_SCRIPTING_DATATABLE_H_

#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>
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
/// Two of them need it. A compiled kernel called once per element must not pay a
/// std::function per point -- measured at 0.4 us per face just to BUILD the
/// table, against 0.3 us to evaluate it, so the indirection already costs more
/// than the arithmetic. And a device kernel cannot call a std::function at all,
/// at any price.
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
  /// Points per element: point p reads element p / divisor -- 1 for a column per point, the
  /// points per cell for a column per cell (bindCellView).
  std::size_t divisor{1};
  /// Optional: element of the point set -> element in memory, for a gathered subset.
  const std::uint32_t* index{nullptr};
  /// One value for all points (bindConstant). A device kernel takes it by value, read on the host
  /// when a call is packed -- from the base of the call if the call moves it -- so its base may be
  /// host memory, and a value that changes from call to call costs no copy to the device.
  bool uniform{false};

  /// The element point `point` reads; base + element * byteStride + byteOffset is its address.
  [[nodiscard]] std::size_t element(std::size_t point) const {
    const std::size_t local = point / divisor;
    return index == nullptr ? local : static_cast<std::size_t>(index[local]);
  }
  /// Whether the view maps points onto elements one to one, in order.
  [[nodiscard]] bool pointwise() const { return divisor == 1 && index == nullptr; }
};

/// A vector of coefficients per cell, for one named quantity -- what a contraction
/// (expr::Kind::Contract) reads instead of a column. Coefficient m of cell c lies at
///   base + index(c) * cellStride + m * modeStride   (in bytes),
/// with index(c) = cellIndex[c] if given and c otherwise; the cells are those of the point set,
/// point p belonging to cell p / pointsPerCell.
///
/// One input per quantity: an interleaved layout dofs[cell][mode][quantity] binds quantity q at
/// base + q * paddedModes * sizeof(T) with modeStride = sizeof(T), and a fused simulation s at
/// base + s * sizeof(T) with modeStride = NumSimulations * sizeof(T).
struct BlockInput {
  const void* base{nullptr};
  std::size_t cellStride{0};
  std::size_t modeStride{0};
  /// Cell of the point set -> cell in memory, for a gathered subset of cells; null is the identity.
  const std::uint32_t* cellIndex{nullptr};
  DataType type{DataType::F64};
  /// Coefficients available per cell; at least as many as a program contracts over.
  std::size_t length{0};
};

/// A matrix blocks are contracted against: entry (row, m) at base[row + m * leadingDimension], one
/// row per point of a cell and one column per coefficient -- the layout of
/// numerical::projection::Table.
struct MatrixInput {
  const void* base{nullptr};
  DataType type{DataType::F64};
  std::size_t rows{0};
  std::size_t cols{0};
  std::size_t leadingDimension{0};
};

/// Where a declared state of a program lives when the consumer keeps it: per cell of the point
/// set, like a block. The state of point p lies at
///   base + index(c) * cellStride + (p - c * pointsPerCell) * pointStride   (in bytes),
/// with c = p / pointsPerCell and index(c) = cellIndex[c] if given and c otherwise, stored in the
/// compute type of the program. A state that no table keeps lives in the Binding instead.
struct StateInput {
  void* base{nullptr};
  std::size_t pointsPerCell{1};
  std::size_t cellStride{0};
  std::size_t pointStride{0};
  /// Cell of the point set -> cell in memory, for a gathered subset of cells; null is the identity.
  const std::uint32_t* cellIndex{nullptr};
  DataType type{DataType::F64};

  /// The offset of the state of point `point` from the base, in bytes.
  [[nodiscard]] std::size_t offset(std::size_t point) const {
    const std::size_t cell = point / pointsPerCell;
    const std::size_t memoryCell = cellIndex == nullptr ? cell : cellIndex[cell];
    return memoryCell * cellStride + (point - cell * pointsPerCell) * pointStride;
  }
};

struct BlockEntry {
  std::string name;
  BlockInput input;
};

struct StateEntry {
  std::string name;
  StateInput input;
};

struct MatrixEntry {
  std::string name;
  MatrixInput input;
};

struct DataEntry {
  std::string name;
  Direction direction;
  DataType datatype{DataType::F64};
  std::function<void(std::size_t, void*)> accessor;
  std::function<void(std::size_t, const void*)> setter;
  /// Empty for bindComputed and bindComputedBatch; set by every other bind form.
  std::optional<StridedView> view;
  /// Fills the contiguous range [first, first + count) at once, in the column's own type. Set by
  /// bindComputedBatch only; the per-point `accessor` is then the same function with count 1.
  std::function<void(std::size_t first, std::size_t count, void* out)> batchAccessor;

  template <typename T>
  [[nodiscard]] T getValue(std::size_t index) const {
    assert(DataTypeTraits<T>::Type == datatype);

    T out{};
    accessor(index, &out);
    return out;
  }

  /// The values of [first, first + count), in the column's own type: one call for a batch-computed
  /// column, a strided copy for a view, and one accessor call per point otherwise.
  template <typename T>
  void getValues(std::size_t first, std::size_t count, T* out) const {
    assert(DataTypeTraits<T>::Type == datatype);

    if (batchAccessor) {
      batchAccessor(first, count, out);
    } else if (view.has_value()) {
      const auto* bytes = static_cast<const char*>(view->base) + view->byteOffset;
      for (std::size_t i = 0; i < count; ++i) {
        std::memcpy(out + i, bytes + view->element(first + i) * view->byteStride, sizeof(T));
      }
    } else {
      for (std::size_t i = 0; i < count; ++i) {
        accessor(first + i, out + i);
      }
    }
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

/// A number of its own for every object, from a process-wide counter: a copy gets a new one, and
/// so does an object moved from.
class InstanceId {
  public:
  InstanceId() : value_(next()) {}
  InstanceId(const InstanceId& /*other*/) : value_(next()) {}
  InstanceId(InstanceId&& other) noexcept : value_(other.value_) { other.value_ = next(); }
  InstanceId& operator=(const InstanceId& other) {
    if (this != &other) {
      value_ = next();
    }
    return *this;
  }
  InstanceId& operator=(InstanceId&& other) noexcept {
    if (this != &other) {
      value_ = other.value_;
      other.value_ = next();
    }
    return *this;
  }
  ~InstanceId() = default;

  [[nodiscard]] std::uint64_t value() const { return value_; }

  private:
  static std::uint64_t next() noexcept {
    static std::atomic<std::uint64_t> counter{0};
    return ++counter;
  }

  std::uint64_t value_;
};

class DataTable {
  public:
  explicit DataTable(std::size_t numPoints) : numPoints_(numPoints) {}

  /// The table as it is bound: the same as long as nothing is bound to it, and different for any
  /// other table -- also for a copy, and for one that takes the place of a destroyed one at its
  /// address. A binding resolves against the columns of a table as they were bound, so one that
  /// is kept for later calls is kept for this revision.
  [[nodiscard]] std::pair<std::uint64_t, std::size_t> revision() const {
    return {instance_.value(),
            dataEntries_.size() + blockEntries_.size() + matrixEntries_.size() +
                stateEntries_.size()};
  }

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
                                        makeView(base, stride, offset, true),
                                        nullptr});
  }

  /// A column that is constant per cell of `pointsPerCell` points: point p reads
  /// base[c * stride + offset] with c = p / pointsPerCell, or c = cellIndex[p / pointsPerCell]
  /// for a gathered subset of cells. Input only.
  template <typename T>
  void bindCellView(std::string name,
                    const T* base,
                    std::size_t pointsPerCell,
                    std::size_t stride = 1,
                    std::size_t offset = 0,
                    const std::uint32_t* cellIndex = nullptr) {
    auto view = makeView(base, stride, offset, false);
    view.divisor = pointsPerCell;
    view.index = cellIndex;
    const auto accessor = [=](std::size_t idx, void* out) {
      *reinterpret_cast<T*>(out) = base[view.element(idx) * stride + offset];
    };
    dataEntries_.emplace_back(DataEntry{
        std::move(name), Direction::In, DataTypeTraits<T>::Type, accessor, nullptr, view, nullptr});
  }

  /// A value that is the same at every point, copied when bound: a stride-0 view
  /// onto a copy the table holds, so that a table that merely names a constant
  /// (such as `sim`) can still reach a device kernel, which bindComputed, with
  /// no address arithmetic behind it, cannot. For a value that changes from
  /// call to call, pass its address as the base of the call (KernelArgs::inputs):
  /// the view is uniform, so a device kernel takes the value by value.
  template <typename T>
  void bindConstant(std::string name, const T& value) {
    auto held = std::make_shared<T>(value);
    const auto accessor = [held](std::size_t, void* out) { *reinterpret_cast<T*>(out) = *held; };
    StridedView view;
    view.base = const_cast<void*>(static_cast<const void*>(held.get()));
    view.byteStride = 0;
    view.byteOffset = 0;
    view.writable = false;
    view.uniform = true;
    dataEntries_.emplace_back(DataEntry{
        std::move(name), Direction::In, DataTypeTraits<T>::Type, accessor, nullptr, view, nullptr});
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
                                        makeView(base, stride, offset, false),
                                        nullptr});
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
                                        makeMemberView(base, member, true),
                                        nullptr});
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
                                        makeMemberView(base, member, false),
                                        nullptr});
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
                                        std::nullopt,
                                        nullptr});
  }

  /// Computed, a contiguous range at a time: `fn(first, count, out)` writes the values of the
  /// points [first, first + count) to out[0 .. count). Consumers that work in tiles call it once
  /// per tile, which amortises the call -- and lets the producer share work between neighbouring
  /// points, such as the transformation of a cell across its quadrature points -- without
  /// materialising the whole column. Like bindComputed, it has no address arithmetic behind it and
  /// hence cannot reach a device kernel. `fn` has to be safe to call concurrently on disjoint
  /// ranges.
  template <typename T, typename F>
  void bindComputedBatch(std::string name, F&& fn) {
    auto batch = std::make_shared<std::decay_t<F>>(std::forward<F>(fn));

    const auto accessor = [batch](std::size_t idx, void* out) {
      (*batch)(idx, std::size_t{1}, reinterpret_cast<T*>(out));
    };
    const auto batchAccessor = [batch](std::size_t first, std::size_t count, void* out) {
      (*batch)(first, count, reinterpret_cast<T*>(out));
    };

    dataEntries_.emplace_back(DataEntry{std::move(name),
                                        Direction::In,
                                        DataTypeTraits<T>::Type,
                                        accessor,
                                        nullptr,
                                        std::nullopt,
                                        batchAccessor});
  }

  /// A block of `length` coefficients per cell; strides in elements of T (cf. BlockInput).
  template <typename T>
  void bindBlock(std::string name,
                 const T* base,
                 std::size_t length,
                 std::size_t cellStride,
                 std::size_t modeStride = 1,
                 const std::uint32_t* cellIndex = nullptr) {
    BlockInput input;
    input.base = base;
    input.cellStride = cellStride * sizeof(T);
    input.modeStride = modeStride * sizeof(T);
    input.cellIndex = cellIndex;
    input.type = DataTypeTraits<T>::Type;
    input.length = length;
    blockEntries_.push_back(BlockEntry{std::move(name), input});
  }

  void bindBlock(std::string name, const BlockInput& input) {
    blockEntries_.push_back(BlockEntry{std::move(name), input});
  }

  /// A matrix with `rows` points per cell and `cols` coefficients (cf. MatrixInput).
  template <typename T>
  void bindMatrix(std::string name,
                  const T* base,
                  std::size_t rows,
                  std::size_t cols,
                  std::size_t leadingDimension) {
    MatrixInput input;
    input.base = base;
    input.type = DataTypeTraits<T>::Type;
    input.rows = rows;
    input.cols = cols;
    input.leadingDimension = leadingDimension;
    matrixEntries_.push_back(MatrixEntry{std::move(name), input});
  }

  [[nodiscard]] std::size_t numPoints() const { return numPoints_; }

  [[nodiscard]] const std::vector<DataEntry>& dataEntries() const { return dataEntries_; }
  /// Keeps the state `name` of a program at `base` (cf. StateInput); strides in elements of T.
  template <typename T>
  void bindState(std::string name,
                 T* base,
                 std::size_t pointsPerCell,
                 std::size_t cellStride,
                 std::size_t pointStride,
                 const std::uint32_t* cellIndex = nullptr) {
    StateInput input;
    input.base = base;
    input.pointsPerCell = pointsPerCell;
    input.cellStride = cellStride * sizeof(T);
    input.pointStride = pointStride * sizeof(T);
    input.cellIndex = cellIndex;
    input.type = DataTypeTraits<T>::Type;
    stateEntries_.push_back(StateEntry{std::move(name), input});
  }

  void bindState(std::string name, const StateInput& input) {
    stateEntries_.push_back(StateEntry{std::move(name), input});
  }

  [[nodiscard]] const std::vector<BlockEntry>& blockEntries() const { return blockEntries_; }
  [[nodiscard]] const std::vector<StateEntry>& stateEntries() const { return stateEntries_; }
  [[nodiscard]] const std::vector<MatrixEntry>& matrixEntries() const { return matrixEntries_; }

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
  std::vector<BlockEntry> blockEntries_;
  std::vector<MatrixEntry> matrixEntries_;
  std::vector<StateEntry> stateEntries_;
  InstanceId instance_;
};

} // namespace seissol::reader::scripting
#endif // SEISSOL_SRC_READER_SCRIPTING_DATATABLE_H_
