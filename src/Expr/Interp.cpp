// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/Interp.h"

#include "Expr/Expand.h"
#include "Expr/Ir.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Interpolation.h"
#include "Reader/Scripting/DataTable.h"
#include "utils/logger.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace seissol::expr {

namespace {

[[noreturn]] void fail(const std::string& what) { throw std::invalid_argument("expr: " + what); }

using reader::scripting::DataType;

template <typename T, typename MatrixT, typename BlockT>
void contractTyped(const ContractOperands& operands,
                   const std::size_t* points,
                   std::size_t first,
                   std::size_t count,
                   T* dst) {
  const auto* matrix = static_cast<const MatrixT*>(operands.matrix);
  const auto* block = static_cast<const char*>(operands.block);
  const std::size_t rows = operands.rows;
  const std::size_t cols = operands.cols;
  const std::size_t leadingDimension = operands.leadingDimension;
  for (std::size_t lane = 0; lane < count; ++lane) {
    const std::size_t point = points == nullptr ? first + lane : points[lane];
    const std::size_t cell = point / rows;
    const std::size_t row = point - cell * rows;
    const std::size_t memoryCell =
        operands.cellIndex == nullptr ? cell : static_cast<std::size_t>(operands.cellIndex[cell]);
    const char* coefficients = block + memoryCell * operands.cellStride;
    const MatrixT* matrixRow = matrix + row;
    T sum = T(0);
    for (std::size_t m = 0; m < cols; ++m) {
      BlockT coefficient{};
      std::memcpy(&coefficient, coefficients + m * operands.modeStride, sizeof(BlockT));
      sum += static_cast<T>(matrixRow[m * leadingDimension]) * static_cast<T>(coefficient);
    }
    dst[lane] = sum;
  }
}

template <typename T, typename MatrixT>
void contractBlockTyped(const ContractOperands& operands,
                        const std::size_t* points,
                        std::size_t first,
                        std::size_t count,
                        T* dst) {
  switch (operands.blockType) {
  case DataType::F32:
    contractTyped<T, MatrixT, float>(operands, points, first, count, dst);
    return;
  case DataType::F64:
    contractTyped<T, MatrixT, double>(operands, points, first, count, dst);
    return;
  default:
    fail("a contraction block has to be stored as f32 or f64");
  }
}

} // namespace

template <typename T>
void contractLanes(const ContractOperands& operands,
                   const std::size_t* points,
                   std::size_t first,
                   std::size_t count,
                   T* dst) {
  switch (operands.matrixType) {
  case DataType::F32:
    contractBlockTyped<T, float>(operands, points, first, count, dst);
    return;
  case DataType::F64:
    contractBlockTyped<T, double>(operands, points, first, count, dst);
    return;
  default:
    fail("a contraction matrix has to be stored as f32 or f64");
  }
}

template void contractLanes<double>(
    const ContractOperands&, const std::size_t*, std::size_t, std::size_t, double*);
template void contractLanes<float>(
    const ContractOperands&, const std::size_t*, std::size_t, std::size_t, float*);

template <typename T>
ContractOperands TileIo<T>::contraction(MatrixId /*matrix*/, BlockId /*block*/) const {
  fail("the program contracts a block, but its tile source binds no blocks");
}

template <typename T>
void TileIo<T>::gatherState(std::size_t /*first*/, std::size_t /*count*/, T* /*dst*/) const {
  fail("the program declares state, but its tile source keeps none");
}

template <typename T>
void TileIo<T>::scatterState(std::size_t /*first*/, std::size_t /*count*/, const T* /*src*/) {
  fail("the program declares state, but its tile source keeps none");
}

std::size_t chooseTileSize(std::int32_t peakSlots, ComputeType type, std::size_t budgetBytes) {
  const std::size_t element = (type == ComputeType::F32) ? sizeof(float) : sizeof(double);
  const auto slots = static_cast<std::size_t>(std::max<std::int32_t>(peakSlots, 1));
  std::size_t lanes = budgetBytes / (slots * element);
  lanes = (lanes / TileLaneGranularity) * TileLaneGranularity;
  return std::clamp(lanes, MinTileSize, MaxTileSize);
}

namespace {

/// The bytes of a transient of `type` per lane: a truth value takes one (0 or 1).
std::size_t widthOf(ComputeType type) {
  switch (type) {
  case ComputeType::F32:
    return sizeof(float);
  case ComputeType::F64:
    return sizeof(double);
  case ComputeType::I64:
    return sizeof(std::int64_t);
  case ComputeType::Bool:
    return sizeof(std::uint8_t);
  }
  return sizeof(double);
}

/// The offsets of the transients of `stage` in the scratch of a tile of `lanes` lanes, in bytes,
/// each slot on a cache line of its own, whatever the type of the one before; and their end.
std::size_t slotOffsets(const StageCode& stage, std::size_t lanes, std::vector<std::size_t>& out) {
  constexpr std::size_t Line = 64;
  out.clear();
  std::size_t offset = 0;
  for (std::int32_t slot = 0; slot < stage.slotCount; ++slot) {
    out.push_back(offset);
    const auto type = static_cast<std::size_t>(slot) < stage.slotTypes.size()
                          ? stage.slotTypes[static_cast<std::size_t>(slot)]
                          : ComputeType::F64;
    offset += (widthOf(type) * lanes + Line - 1) / Line * Line;
  }
  return offset;
}

/// The transients of the larger stage in slots of the compute type: what the tile size follows.
std::int32_t peakSlotsIn(const LoweredProgram& lowered, ComputeType computeType) {
  std::size_t bytes = 0;
  for (const StageCode* stage : {&lowered.precompute(), &lowered.run()}) {
    std::size_t perLane = 0;
    for (std::int32_t slot = 0; slot < stage->slotCount; ++slot) {
      const auto type = static_cast<std::size_t>(slot) < stage->slotTypes.size()
                            ? stage->slotTypes[static_cast<std::size_t>(slot)]
                            : ComputeType::F64;
      perLane += widthOf(type);
    }
    bytes = std::max(bytes, perLane);
  }
  const std::size_t width = widthOf(computeType);
  return static_cast<std::int32_t>((bytes + width - 1) / width);
}

/// The lanes of a pointwise operation in the type U.
template <typename U>
void pwLanes(Fn fn,
             const U* __restrict a0,
             const U* __restrict a1,
             const U* __restrict a2,
             U* __restrict dst,
             std::size_t count) {
  using T = U;
  switch (fn) {
#define SEISSOL_EXPR_PW_LOOP_UNARY(NAME, EXPR)                                                     \
  case Fn::NAME: {                                                                                 \
    _Pragma("omp simd") for (std::size_t l = 0; l < count; ++l) {                                  \
      const T x = a0[l];                                                                           \
      dst[l] = (EXPR);                                                                             \
    }                                                                                              \
    break;                                                                                         \
  }
#define SEISSOL_EXPR_PW_LOOP_BINARY(NAME, EXPR)                                                    \
  case Fn::NAME: {                                                                                 \
    _Pragma("omp simd") for (std::size_t l = 0; l < count; ++l) {                                  \
      const T x = a0[l];                                                                           \
      const T y = a1[l];                                                                           \
      dst[l] = (EXPR);                                                                             \
    }                                                                                              \
    break;                                                                                         \
  }
#define SEISSOL_EXPR_PW_LOOP_TERNARY(NAME, EXPR)                                                   \
  case Fn::NAME: {                                                                                 \
    _Pragma("omp simd") for (std::size_t l = 0; l < count; ++l) {                                  \
      const T x = a0[l];                                                                           \
      const T y = a1[l];                                                                           \
      const T z = a2[l];                                                                           \
      dst[l] = (EXPR);                                                                             \
    }                                                                                              \
    break;                                                                                         \
  }
    SEISSOL_EXPR_PW_LIST(
        SEISSOL_EXPR_PW_LOOP_UNARY, SEISSOL_EXPR_PW_LOOP_BINARY, SEISSOL_EXPR_PW_LOOP_TERNARY)
#undef SEISSOL_EXPR_PW_LOOP_UNARY
#undef SEISSOL_EXPR_PW_LOOP_BINARY
#undef SEISSOL_EXPR_PW_LOOP_TERNARY
  }
}

/// The lanes of a pointwise operation on whole numbers: those the written-out lookups compute (the
/// minimum compares whole numbers, where the table's fmin would compute in binary64).
void wholeLanes(Fn fn,
                const std::int64_t* __restrict a0,
                const std::int64_t* __restrict a1,
                std::int64_t* __restrict dst,
                std::size_t count) {
  switch (fn) {
  case Fn::Neg:
#pragma omp simd
    for (std::size_t l = 0; l < count; ++l) {
      dst[l] = -a0[l];
    }
    return;
  case Fn::Add:
#pragma omp simd
    for (std::size_t l = 0; l < count; ++l) {
      dst[l] = a0[l] + a1[l];
    }
    return;
  case Fn::Sub:
#pragma omp simd
    for (std::size_t l = 0; l < count; ++l) {
      dst[l] = a0[l] - a1[l];
    }
    return;
  case Fn::Mul:
#pragma omp simd
    for (std::size_t l = 0; l < count; ++l) {
      dst[l] = a0[l] * a1[l];
    }
    return;
  case Fn::Min:
#pragma omp simd
    for (std::size_t l = 0; l < count; ++l) {
      dst[l] = a1[l] < a0[l] ? a1[l] : a0[l];
    }
    return;
  case Fn::Max:
#pragma omp simd
    for (std::size_t l = 0; l < count; ++l) {
      dst[l] = a0[l] < a1[l] ? a1[l] : a0[l];
    }
    return;
  default:
    fail(std::string("`") + name(fn) + "` does not compute on whole numbers");
  }
}

/// dst[l] = U(src[l]) for the lanes, from data of type V; into a truth value, whether it is not 0.
template <typename U, typename V>
void convertLanes(const V* __restrict src, U* __restrict dst, std::size_t count) {
  if constexpr (std::is_same_v<U, std::uint8_t>) {
#pragma omp simd
    for (std::size_t l = 0; l < count; ++l) {
      dst[l] = src[l] != V(0) ? 1 : 0;
    }
  } else {
#pragma omp simd
    for (std::size_t l = 0; l < count; ++l) {
      dst[l] = static_cast<U>(src[l]);
    }
  }
}

/// dst[l] = the element of the array `data` at offsets[l] + shift.
template <typename U, typename E>
void gatherLanes(const void* data,
                 const std::int64_t* __restrict offsets,
                 std::int64_t shift,
                 U* __restrict dst,
                 std::size_t count) {
  const auto* elements = static_cast<const E*>(data);
  for (std::size_t l = 0; l < count; ++l) {
    dst[l] = static_cast<U>(static_cast<double>(elements[offsets[l] + shift]));
  }
}

/// The truth of `fn` (Lt, Le, Eq) for the lanes of two operands of type V.
template <typename V>
void compareLanes(Fn fn,
                  const V* __restrict a,
                  const V* __restrict b,
                  std::uint8_t* __restrict dst,
                  std::size_t count) {
  switch (fn) {
  case Fn::Lt:
#pragma omp simd
    for (std::size_t l = 0; l < count; ++l) {
      dst[l] = a[l] < b[l] ? 1 : 0;
    }
    return;
  case Fn::Le:
#pragma omp simd
    for (std::size_t l = 0; l < count; ++l) {
      dst[l] = a[l] <= b[l] ? 1 : 0;
    }
    return;
  case Fn::Eq:
#pragma omp simd
    for (std::size_t l = 0; l < count; ++l) {
      dst[l] = a[l] == b[l] ? 1 : 0;
    }
    return;
  default:
    fail(std::string("`") + name(fn) + "` gives no truth value");
  }
}

} // namespace

template <typename T>
TileInterpreter<T>::TileInterpreter(const Program& program,
                                    const LoweredProgram& lowered,
                                    const GridArrays* grids,
                                    const InterpreterOptions& options)
    : lowered_(&lowered), grids_(grids), perturbation_(options.perturbation),
      decisions_(options.decisions), branches_(options.branches),
      numInputs_(static_cast<std::int32_t>(program.inputs().size())),
      numOutputs_(static_cast<std::int32_t>(program.outputs().size())),
      stateSlots_(lowered.stateSlotCount()) {
  // the data moves in its type; the operations compute in theirs
  const bool wantsF32 = lowered.dataType() == ComputeType::F32;
  if (wantsF32 != std::is_same_v<T, float>) {
    fail("interpreter instantiated for another type than the program moves its data in");
  }
  tileSize_ = options.tileSize != 0 ? options.tileSize
                                    : chooseTileSize(peakSlotsIn(lowered, program.computeType()),
                                                     program.computeType(),
                                                     options.budgetBytes);
  const std::size_t precomputeBytes =
      slotOffsets(lowered.precompute(), tileSize_, precomputeOffsets_);
  const std::size_t runBytes = slotOffsets(lowered.run(), tileSize_, runOffsets_);
  scratch_.assign(std::max<std::size_t>(std::max(precomputeBytes, runBytes), sizeof(double)) /
                          sizeof(double) +
                      1,
                  0.0);
  inputTile_.assign(static_cast<std::size_t>(std::max(numInputs_, 1)) * tileSize_, T(0));
  outputTile_.assign(static_cast<std::size_t>(std::max(numOutputs_, 1)) * tileSize_, T(0));
  stateTile_.assign(static_cast<std::size_t>(std::max(stateSlots_, 1)) * tileSize_, T(0));
}

// One switch per instruction, amortised over `count` lanes. The lane loops are
// left in the plainest possible shape — unit stride, no aliasing, no branches
// except the ones inside a Select's arithmetic — so the vectorizer can take
// them. Each instruction computes in its type, its operands in theirs; the data
// arrives and leaves in T.
template <typename T>
void TileInterpreter<T>::runStage(const StageCode& stage,
                                  const TileIo<T>& io,
                                  const T* inputTile,
                                  T* outputTile,
                                  T* stateTile,
                                  T* persistent,
                                  std::size_t numPoints,
                                  std::size_t first,
                                  std::size_t count) {
  const std::vector<std::int32_t>& operands = lowered_->operands();
  auto* const scratchBase = reinterpret_cast<unsigned char*>(scratch_.data());
  const std::vector<std::size_t>& offsets =
      &stage == &lowered_->precompute() ? precomputeOffsets_ : runOffsets_;
  const auto typeOf = [&](std::int32_t slot) {
    return static_cast<std::size_t>(slot) < stage.slotTypes.size()
               ? stage.slotTypes[static_cast<std::size_t>(slot)]
               : ComputeType::F64;
  };
  const auto bytesOf = [&](std::int32_t slot) {
    return scratchBase + offsets[static_cast<std::size_t>(slot)];
  };

  // the lanes of the slot `slot`, in its type, converted into `dst`
  const auto convertSlot = [&](std::int32_t slot, auto* dst) {
    switch (typeOf(slot)) {
    case ComputeType::F32:
      convertLanes(reinterpret_cast<const float*>(bytesOf(slot)), dst, count);
      return;
    case ComputeType::F64:
      convertLanes(reinterpret_cast<const double*>(bytesOf(slot)), dst, count);
      return;
    case ComputeType::I64:
      convertLanes(reinterpret_cast<const std::int64_t*>(bytesOf(slot)), dst, count);
      return;
    case ComputeType::Bool:
      convertLanes(reinterpret_cast<const std::uint8_t*>(bytesOf(slot)), dst, count);
      return;
    }
  };
  // the slot `slot` into the data of type T at `dst`
  const auto store = [&](std::int32_t slot, T* dst) { convertSlot(slot, dst); };

  const auto arrayOf = [&](GridId grid) -> const reader::datafield::ArrayView& {
    if (grids_ == nullptr) {
      fail("program samples a grid, but no grid arrays were supplied");
    }
    return grids_->array(grid);
  };

  // the signature of the decisions of the selects of each point, where asked for
  const auto decide = [&](const Instruction& inst, const auto& taken) {
    const auto index = static_cast<std::uint64_t>(&inst - stage.code.data());
    if (decisions_ != nullptr) {
      for (std::size_t l = 0; l < count; ++l) {
        auto& signature = (*decisions_)[first + l];
        signature = splitMix64(signature ^ ((index << 1U) | (taken(l) ? 1U : 0U)));
      }
    }
    if (branches_ != nullptr && &stage == &lowered_->run()) {
      auto& branch = (*branches_)[static_cast<std::size_t>(index)];
      branch.resize(numPoints, 0);
      for (std::size_t l = 0; l < count; ++l) {
        branch[first + l] = taken(l) ? 1 : 0;
      }
    }
  };

  // a product, or an operation of the maths library, moved by its ulps up or down per point
  const auto perturb = [&](const Instruction& inst, auto* dst) {
    using U = std::remove_pointer_t<decltype(dst)>;
    int ulps = 0;
    switch (inst.fn) {
    case Fn::Mul:
      ulps = inst.strict ? 0 : perturbation_.productUlps;
      break;
    case Fn::Exp:
    case Fn::Log:
    case Fn::Log2:
    case Fn::Log10:
    case Fn::Sin:
    case Fn::Cos:
    case Fn::Tan:
    case Fn::Asin:
    case Fn::Acos:
    case Fn::Atan:
    case Fn::Atan2:
    case Fn::Sinh:
    case Fn::Cosh:
    case Fn::Tanh:
    case Fn::Asinh:
    case Fn::Acosh:
    case Fn::Atanh:
    case Fn::Erf:
    case Fn::Pow:
      ulps = perturbation_.boundedUlps;
      break;
    default:
      break;
    }
    if (ulps == 0) {
      return;
    }
    const auto index = static_cast<std::uint64_t>(&inst - stage.code.data());
    for (std::size_t l = 0; l < count; ++l) {
      const std::uint64_t draw =
          splitMix64(perturbation_.seed ^ splitMix64(((first + l) << 20U) ^ index));
      const U toward = (draw & 1U) != 0 ? std::numeric_limits<U>::infinity()
                                        : -std::numeric_limits<U>::infinity();
      for (int step = 0; step < ulps; ++step) {
        dst[l] = std::nextafter(dst[l], toward);
      }
    }
  };

  const auto execute = [&](const Instruction& inst, auto tag) {
    using U = decltype(tag);
    U* const __restrict dst = reinterpret_cast<U*>(bytesOf(inst.dst));
    const auto operand = [&](std::int32_t k) {
      return reinterpret_cast<const U*>(bytesOf(operands[inst.operandBegin + k]));
    };
    switch (inst.op) {
    case Opcode::Const: {
      const U value = static_cast<U>(inst.value);
#pragma omp simd
      for (std::size_t l = 0; l < count; ++l) {
        dst[l] = value;
      }
      break;
    }
    case Opcode::LoadInput:
      convertLanes(inputTile + static_cast<std::size_t>(inst.imm) * count, dst, count);
      break;
    case Opcode::LoadPersistent:
      // a declared state from the tile, a hoisted value from the persistent buffer
      convertLanes(inst.imm < stateSlots_
                       ? stateTile + static_cast<std::size_t>(inst.imm) * count
                       : persistent + static_cast<std::size_t>(inst.imm - stateSlots_) * numPoints +
                             first,
                   dst,
                   count);
      break;
    case Opcode::GridParam: {
      // read when the stage runs: the window of a time-dependent grid moves it between calls
      const U value = static_cast<U>(
          gridField(arrayOf(inst.grid), static_cast<GridField>(inst.comp), inst.imm));
#pragma omp simd
      for (std::size_t l = 0; l < count; ++l) {
        dst[l] = value;
      }
      break;
    }
    case Opcode::GridLoad: {
      const reader::datafield::ArrayView& view = arrayOf(inst.grid);
      const std::int32_t source = operands[inst.operandBegin];
      if (typeOf(source) != ComputeType::I64) {
        fail("the offset of an element of a grid is no whole number");
      }
      const auto* offsets = reinterpret_cast<const std::int64_t*>(bytesOf(source));
      if (inst.comp < 0) {
        // a knot, or the reciprocal of a spacing, of an axis of knots (knotsOf())
        const auto axis = static_cast<std::size_t>((-1 - inst.comp) / 2);
        const bool knots = (-1 - inst.comp) % 2 == 0;
        const double* array = knots ? view.knots[axis] : view.inverseSpacings[axis];
        if (array == nullptr) {
          fail("an axis of a grid has no knots");
        }
        gatherLanes<U, double>(array, offsets, 0, dst, count);
        break;
      }
      const auto shift =
          static_cast<std::int64_t>(inst.comp) * static_cast<std::int64_t>(view.componentStride);
      if (view.type == reader::datafield::ElementType::Float) {
        gatherLanes<U, float>(view.data, offsets, shift, dst, count);
      } else {
        gatherLanes<U, double>(view.data, offsets, shift, dst, count);
      }
      break;
    }
    case Opcode::Contract:
      if constexpr (std::is_integral_v<U>) {
        fail("a contraction computes on whole numbers");
      } else {
        contractLanes<U>(
            io.contraction(inst.matrix, inst.block), io.points(first), first, count, dst);
      }
      break;
    case Opcode::Pw: {
      if (inst.fn == Fn::Select && typeOf(operands[inst.operandBegin]) == ComputeType::Bool) {
        // a truth value chooses between operands of the type of the instruction
        const auto* condition =
            reinterpret_cast<const std::uint8_t*>(bytesOf(operands[inst.operandBegin]));
        const U* const a1 = operand(1);
        const U* const a2 = operand(2);
#pragma omp simd
        for (std::size_t l = 0; l < count; ++l) {
          dst[l] = condition[l] != 0 ? a1[l] : a2[l];
        }
        decide(inst, [&](std::size_t l) { return condition[l] != 0; });
        break;
      }
      const U* const a0 = operand(0);
      const U* const a1 = inst.operandCount > 1 ? operand(1) : a0;
      const U* const a2 = inst.operandCount > 2 ? operand(2) : a0;
      if constexpr (std::is_integral_v<U>) {
        wholeLanes(inst.fn, a0, a1, dst, count);
      } else {
        pwLanes<U>(inst.fn, a0, a1, a2, dst, count);
        if (inst.fn == Fn::Select) {
          decide(inst, [&](std::size_t l) { return a0[l] != U(0); });
        }
        if (perturbation_.on) {
          perturb(inst, dst);
        }
      }
      break;
    }
    case Opcode::Convert:
      convertSlot(operands[inst.operandBegin], dst);
      break;
    }
  };

  // An instruction of type Bool: a comparison of operands in their own type, `and` and `or` of
  // truth values, or a truth value chosen, converted, loaded or constant.
  const auto executeBool = [&](const Instruction& inst) {
    auto* const __restrict dst = reinterpret_cast<std::uint8_t*>(bytesOf(inst.dst));
    const auto lanes = [&](std::int32_t k) { return bytesOf(operands[inst.operandBegin + k]); };
    switch (inst.op) {
    case Opcode::Const: {
      const std::uint8_t value = inst.value != 0.0 ? 1 : 0;
      for (std::size_t l = 0; l < count; ++l) {
        dst[l] = value;
      }
      return;
    }
    case Opcode::LoadPersistent:
      convertLanes(inst.imm < stateSlots_
                       ? stateTile + static_cast<std::size_t>(inst.imm) * count
                       : persistent + static_cast<std::size_t>(inst.imm - stateSlots_) * numPoints +
                             first,
                   dst,
                   count);
      return;
    case Opcode::Convert:
      convertSlot(operands[inst.operandBegin], dst);
      return;
    case Opcode::Pw:
      switch (inst.fn) {
      case Fn::Lt:
      case Fn::Le:
      case Fn::Eq:
        switch (typeOf(operands[inst.operandBegin])) {
        case ComputeType::F32:
          compareLanes(inst.fn,
                       reinterpret_cast<const float*>(lanes(0)),
                       reinterpret_cast<const float*>(lanes(1)),
                       dst,
                       count);
          return;
        case ComputeType::F64:
          compareLanes(inst.fn,
                       reinterpret_cast<const double*>(lanes(0)),
                       reinterpret_cast<const double*>(lanes(1)),
                       dst,
                       count);
          return;
        case ComputeType::I64:
          compareLanes(inst.fn,
                       reinterpret_cast<const std::int64_t*>(lanes(0)),
                       reinterpret_cast<const std::int64_t*>(lanes(1)),
                       dst,
                       count);
          return;
        case ComputeType::Bool:
          compareLanes(inst.fn,
                       reinterpret_cast<const std::uint8_t*>(lanes(0)),
                       reinterpret_cast<const std::uint8_t*>(lanes(1)),
                       dst,
                       count);
          return;
        }
        return;
      case Fn::And:
      case Fn::Or: {
        const auto* x = reinterpret_cast<const std::uint8_t*>(lanes(0));
        const auto* y = reinterpret_cast<const std::uint8_t*>(lanes(1));
        if (inst.fn == Fn::And) {
#pragma omp simd
          for (std::size_t l = 0; l < count; ++l) {
            dst[l] = (x[l] != 0 && y[l] != 0) ? 1 : 0;
          }
        } else {
#pragma omp simd
          for (std::size_t l = 0; l < count; ++l) {
            dst[l] = (x[l] != 0 || y[l] != 0) ? 1 : 0;
          }
        }
        return;
      }
      case Fn::Select: {
        const auto* condition = reinterpret_cast<const std::uint8_t*>(lanes(0));
        const auto* x = reinterpret_cast<const std::uint8_t*>(lanes(1));
        const auto* y = reinterpret_cast<const std::uint8_t*>(lanes(2));
#pragma omp simd
        for (std::size_t l = 0; l < count; ++l) {
          dst[l] = condition[l] != 0 ? x[l] : y[l];
        }
        decide(inst, [&](std::size_t l) { return condition[l] != 0; });
        return;
      }
      default:
        break;
      }
      break;
    default:
      break;
    }
    fail(std::string("an instruction of type bool computes `") + name(inst.op) + "`" +
         (inst.op == Opcode::Pw ? std::string(" ") + name(inst.fn) : std::string()));
  };

  for (const Instruction& inst : stage.code) {
    switch (inst.type) {
    case ComputeType::F32:
      execute(inst, 0.0F);
      break;
    case ComputeType::F64:
      execute(inst, 0.0);
      break;
    case ComputeType::I64:
      execute(inst, std::int64_t{0});
      break;
    case ComputeType::Bool:
      executeBool(inst);
      break;
    }
  }

  // A check of a lookup that found what its kind says at a point of the tile: a sample beyond a
  // grid that fails there, or a time the window of a grid did not serve.
  for (const CheckStore& check : stage.checks) {
    if (typeOf(check.source) != ComputeType::Bool) {
      fail("a check of a lookup is no truth value");
    }
    const auto* values = reinterpret_cast<const std::uint8_t*>(bytesOf(check.source));
    std::size_t lane = 0;
    while (lane < count && values[lane] == 0) {
      ++lane;
    }
    if (lane == count) {
      continue;
    }
    if (grids_ == nullptr) {
      fail("a check found a point, but nobody receives it");
    }
    if (check.kind == CheckKind::User) {
      // the first point of the tile that does not meet it, with its inputs
      std::vector<double> inputs(static_cast<std::size_t>(numInputs_));
      for (std::size_t i = 0; i < inputs.size(); ++i) {
        inputs[i] = static_cast<double>(inputTile[i * count + lane]);
      }
      grids_->failed(check.index, first + lane, inputs);
    } else {
      grids_->report(check.grid, check.kind);
    }
  }

  // Stores run after every instruction, which is what makes state updates a
  // parallel assignment rather than a sequential one.
  for (const Store& st : stage.outputs) {
    store(st.source, outputTile + static_cast<std::size_t>(st.target) * count);
  }
  for (const Store& st : stage.persistent) {
    store(st.source,
          st.target < stateSlots_
              ? stateTile + static_cast<std::size_t>(st.target) * count
              : persistent + static_cast<std::size_t>(st.target - stateSlots_) * numPoints + first);
  }
}

namespace {

std::vector<PointRange> resolvePartitions(const std::vector<PointRange>& given,
                                          std::size_t numPoints) {
  if (given.empty()) {
    return {PointRange{0, numPoints}};
  }
  return given;
}

} // namespace

template <typename T>
void TileInterpreter<T>::precompute(const TileIo<T>& io,
                                    std::size_t numPoints,
                                    T* persistent,
                                    const std::vector<PointRange>& partitions) {
  if (!lowered_->hasPrecompute()) {
    return;
  }
  for (const PointRange& range : resolvePartitions(partitions, numPoints)) {
    for (std::size_t first = range.begin; first < range.end; first += tileSize_) {
      const std::size_t count = std::min(tileSize_, range.end - first);
      io.gather(first, count, inputTile_.data());
      // the precompute stage hoists invariant values only, so it never reads a state
      runStage(lowered_->precompute(),
               io,
               inputTile_.data(),
               outputTile_.data(),
               nullptr,
               persistent,
               numPoints,
               first,
               count);
    }
  }
}

template <typename T>
void TileInterpreter<T>::run(TileIo<T>& io,
                             std::size_t numPoints,
                             T* persistent,
                             const std::vector<PointRange>& partitions) {
  for (const PointRange& range : resolvePartitions(partitions, numPoints)) {
    for (std::size_t first = range.begin; first < range.end; first += tileSize_) {
      const std::size_t count = std::min(tileSize_, range.end - first);
      io.gather(first, count, inputTile_.data());
      if (stateSlots_ > 0) {
        io.gatherState(first, count, stateTile_.data());
      }
      runStage(lowered_->run(),
               io,
               inputTile_.data(),
               outputTile_.data(),
               stateTile_.data(),
               persistent,
               numPoints,
               first,
               count);
      io.scatter(first, count, outputTile_.data());
      if (stateSlots_ > 0) {
        io.scatterState(first, count, stateTile_.data());
      }
    }
  }
}

template <typename T>
void initialiseState(const Program& program, T* persistent, std::size_t numPoints) {
  for (std::size_t i = 0; i < program.state().size(); ++i) {
    T* const slot = persistent + i * numPoints;
    const auto value = static_cast<T>(program.state()[i].initial);
    for (std::size_t p = 0; p < numPoints; ++p) {
      slot[p] = value;
    }
  }
}

namespace {

/// The outputs of a tile of one point, and no inputs.
class OutputsOnly final : public TileIo<double> {
  public:
  explicit OutputsOnly(std::vector<double>& outputs) : outputs_(outputs) {}
  void gather(std::size_t /*first*/, std::size_t /*count*/, double* /*dst*/) const override {}
  void scatter(std::size_t /*first*/, std::size_t count, const double* src) override {
    for (std::size_t o = 0; o < outputs_.size(); ++o) {
      outputs_[o] = src[o * count];
    }
  }

  private:
  std::vector<double>& outputs_;
};

} // namespace

std::string interpreterFault() {
  constexpr double Nan = std::numeric_limits<double>::quiet_NaN();
  constexpr double Inf = std::numeric_limits<double>::infinity();
  Program program;
  Arena& arena = program.arena();
  struct Case {
    const char* what;
    NodeId value;
    double expected;
  };
  // 1 + 2^-27 and 1 - 2^-27: their product, 1 - 2^-54, rounds to 1, contracted with - 1 it is not 0
  const NodeId above = arena.konst(1.0 + std::ldexp(1.0, -27));
  const NodeId below = arena.konst(1.0 - std::ldexp(1.0, -27));
  const std::vector<Case> cases = {
      {"a NaN equals itself", arena.pw(Fn::Eq, arena.konst(Nan), arena.konst(Nan)), 0.0},
      {"a NaN is less than 1", arena.pw(Fn::Lt, arena.konst(Nan), arena.konst(1.0)), 0.0},
      {"fmax(NaN, 1) is not 1", arena.pw(Fn::Max, arena.konst(Nan), arena.konst(1.0)), 1.0},
      {"0.3 / 0.1 is computed as 0.3 * (1 / 0.1)",
       arena.pw(Fn::Div, arena.konst(0.3), arena.konst(0.1)),
       std::bit_cast<double>(std::uint64_t{0x4007FFFFFFFFFFFFULL})}, // 2.9999999999999996
      {"a product and a sum are contracted into a fused multiply-add",
       arena.pw(Fn::Add, arena.pw(Fn::Mul, above, below), arena.konst(-1.0)),
       0.0},
      {"an infinity minus itself is not NaN",
       arena.pw(Fn::Sub, arena.konst(Inf), arena.konst(Inf)),
       Nan},
      {"the sign of zero is lost", arena.pw(Fn::Mul, arena.konst(0.0), arena.konst(-1.0)), -0.0}};
  for (std::size_t i = 0; i < cases.size(); ++i) {
    program.addOutput("u" + std::to_string(i), reader::scripting::DataType::F64, cases[i].value);
  }
  const LoweredProgram lowered = lower(program);
  std::vector<double> outputs(cases.size(), 0.5);
  OutputsOnly io(outputs);
  TileInterpreter<double> interpreter(program, lowered, nullptr, InterpreterOptions{});
  std::vector<double> persistent(1, 0.0);
  interpreter.run(io, 1, persistent.data());
  for (std::size_t i = 0; i < cases.size(); ++i) {
    std::uint64_t got = 0;
    std::uint64_t expected = 0;
    std::memcpy(&got, &outputs[i], sizeof(got));
    std::memcpy(&expected, &cases[i].expected, sizeof(expected));
    const bool nan = std::isnan(cases[i].expected);
    if (nan ? !std::isnan(outputs[i]) : got != expected) {
      return cases[i].what;
    }
  }
  return {};
}

void GridArrays::failed(std::int32_t check,
                        std::size_t point,
                        const std::vector<double>& /*inputs*/) const {
  throw std::runtime_error("expr: point " + std::to_string(point) + " does not meet check " +
                           std::to_string(check) + " of the program");
}

template class TileIo<double>;
template class TileIo<float>;
template class TileInterpreter<double>;
template class TileInterpreter<float>;
template void initialiseState<double>(const Program&, double*, std::size_t);
template void initialiseState<float>(const Program&, float*, std::size_t);

} // namespace seissol::expr
