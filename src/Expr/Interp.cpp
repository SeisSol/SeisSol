// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/Interp.h"

#include "Expr/Ir.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Reader/Scripting/DataTable.h"
#include "utils/logger.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
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

/// The bytes of a transient of `type` per lane.
std::size_t widthOf(ComputeType type) {
  return type == ComputeType::F32 ? sizeof(float) : sizeof(double);
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

/// dst[l] = U(src[l]) for the lanes, from data of type V.
template <typename U, typename V>
void convertLanes(const V* __restrict src, U* __restrict dst, std::size_t count) {
#pragma omp simd
  for (std::size_t l = 0; l < count; ++l) {
    dst[l] = static_cast<U>(src[l]);
  }
}

} // namespace

template <typename T>
TileInterpreter<T>::TileInterpreter(const Program& program,
                                    const LoweredProgram& lowered,
                                    const GridSampler* sampler,
                                    const InterpreterOptions& options)
    : lowered_(&lowered), sampler_(sampler),
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

  // the slot `slot` into the data of type T at `dst`
  const auto store = [&](std::int32_t slot, T* dst) {
    if (typeOf(slot) == ComputeType::F32) {
      convertLanes(reinterpret_cast<const float*>(bytesOf(slot)), dst, count);
    } else {
      convertLanes(reinterpret_cast<const double*>(bytesOf(slot)), dst, count);
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
    case Opcode::Lookup: {
      if (sampler_ == nullptr) {
        fail("program contains a grid lookup but no GridSampler was supplied");
      }
      std::array<const U*, MaxLookupDimension> coords{};
      for (std::int32_t k = 0; k < inst.operandCount; ++k) {
        coords[static_cast<std::size_t>(k)] = operand(k);
      }
      sampler_->sampleBatch(inst.grid, inst.comp, coords.data(), inst.operandCount, count, dst);
      break;
    }
    case Opcode::Contract:
      contractLanes<U>(
          io.contraction(inst.matrix, inst.block), io.points(first), first, count, dst);
      break;
    case Opcode::Pw: {
      const U* const a0 = operand(0);
      const U* const a1 = inst.operandCount > 1 ? operand(1) : a0;
      const U* const a2 = inst.operandCount > 2 ? operand(2) : a0;
      pwLanes<U>(inst.fn, a0, a1, a2, dst, count);
      break;
    }
    case Opcode::Convert: {
      const std::int32_t source = operands[inst.operandBegin];
      if (typeOf(source) == ComputeType::F32) {
        convertLanes(reinterpret_cast<const float*>(bytesOf(source)), dst, count);
      } else {
        convertLanes(reinterpret_cast<const double*>(bytesOf(source)), dst, count);
      }
      break;
    }
    }
  };

  for (const Instruction& inst : stage.code) {
    if (inst.type == ComputeType::F32) {
      execute(inst, 0.0F);
    } else {
      execute(inst, 0.0);
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

template class TileIo<double>;
template class TileIo<float>;
template class TileInterpreter<double>;
template class TileInterpreter<float>;
template void initialiseState<double>(const Program&, double*, std::size_t);
template void initialiseState<float>(const Program&, float*, std::size_t);

} // namespace seissol::expr
