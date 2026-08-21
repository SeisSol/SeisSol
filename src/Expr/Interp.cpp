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
#include "utils/logger.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace seissol::expr {

namespace {

[[noreturn]] void fail(const std::string& what) { throw std::invalid_argument("expr: " + what); }

} // namespace

std::size_t chooseTileSize(std::int32_t peakSlots, ComputeType type, std::size_t budgetBytes) {
  const std::size_t element = (type == ComputeType::F32) ? sizeof(float) : sizeof(double);
  const auto slots = static_cast<std::size_t>(std::max<std::int32_t>(peakSlots, 1));
  std::size_t lanes = budgetBytes / (slots * element);
  lanes = (lanes / TileLaneGranularity) * TileLaneGranularity;
  return std::clamp(lanes, MinTileSize, MaxTileSize);
}

template <typename T>
TileInterpreter<T>::TileInterpreter(const Program& program,
                                    const LoweredProgram& lowered,
                                    const GridSampler* sampler,
                                    const InterpreterOptions& options)
    : lowered_(&lowered), sampler_(sampler),
      numInputs_(static_cast<std::int32_t>(program.inputs().size())),
      numOutputs_(static_cast<std::int32_t>(program.outputs().size())) {
  const bool wantsF32 = program.computeType() == ComputeType::F32;
  if (wantsF32 != std::is_same_v<T, float>) {
    fail("interpreter instantiated for the wrong compute type");
  }
  tileSize_ =
      options.tileSize != 0
          ? options.tileSize
          : chooseTileSize(lowered.peakSlotCount(), program.computeType(), options.budgetBytes);
  scratch_.assign(static_cast<std::size_t>(std::max(lowered.peakSlotCount(), 1)) * tileSize_, T(0));
  inputTile_.assign(static_cast<std::size_t>(std::max(numInputs_, 1)) * tileSize_, T(0));
  outputTile_.assign(static_cast<std::size_t>(std::max(numOutputs_, 1)) * tileSize_, T(0));
}

// One switch per instruction, amortised over `count` lanes. The lane loops are
// left in the plainest possible shape — unit stride, no aliasing, no branches
// except the ones inside a Select's arithmetic — so the vectoriser can take
// them. The sderiv interpreter puts `omp simd` on exactly these loops; that is
// still applicable and orthogonal to what happens here.
template <typename T>
void TileInterpreter<T>::runStage(const StageCode& stage,
                                  const T* inputTile,
                                  T* outputTile,
                                  T* persistent,
                                  std::size_t numPoints,
                                  std::size_t first,
                                  std::size_t count) {
  const std::vector<std::int32_t>& operands = lowered_->operands();
  T* const scratchBase = scratch_.data();
  const std::size_t stride = tileSize_;

  const auto slotPtr = [&](std::int32_t slot) { return scratchBase + slot * stride; };

  for (const Instruction& inst : stage.code) {
    T* const __restrict dst = slotPtr(inst.dst);

    switch (inst.op) {
    case Opcode::Const: {
      const T value = static_cast<T>(inst.value);
      for (std::size_t l = 0; l < count; ++l) {
        dst[l] = value;
      }
      break;
    }
    case Opcode::LoadInput: {
      const T* const __restrict src = inputTile + static_cast<std::size_t>(inst.imm) * count;
      for (std::size_t l = 0; l < count; ++l) {
        dst[l] = src[l];
      }
      break;
    }
    case Opcode::LoadPersistent: {
      const T* const __restrict src =
          persistent + static_cast<std::size_t>(inst.imm) * numPoints + first;
      for (std::size_t l = 0; l < count; ++l) {
        dst[l] = src[l];
      }
      break;
    }
    case Opcode::Lookup: {
      if (sampler_ == nullptr) {
        fail("program contains a grid lookup but no GridSampler was supplied");
      }
      std::array<const T*, MaxLookupDimension> coords{};
      for (std::int32_t k = 0; k < inst.operandCount; ++k) {
        coords[static_cast<std::size_t>(k)] = slotPtr(operands[inst.operandBegin + k]);
      }
      sampler_->sampleBatch(inst.grid, inst.comp, coords.data(), inst.operandCount, count, dst);
      break;
    }
    case Opcode::Pw: {
      const T* const __restrict a0 = slotPtr(operands[inst.operandBegin]);
      const T* const __restrict a1 =
          inst.operandCount > 1 ? slotPtr(operands[inst.operandBegin + 1]) : a0;
      const T* const __restrict a2 =
          inst.operandCount > 2 ? slotPtr(operands[inst.operandBegin + 2]) : a0;
      switch (inst.fn) {
#define SEISSOL_EXPR_PW_LOOP_UNARY(NAME, EXPR)                                                     \
  case Fn::NAME: {                                                                                 \
    for (std::size_t l = 0; l < count; ++l) {                                                      \
      const T x = a0[l];                                                                           \
      dst[l] = (EXPR);                                                                             \
    }                                                                                              \
    break;                                                                                         \
  }
#define SEISSOL_EXPR_PW_LOOP_BINARY(NAME, EXPR)                                                    \
  case Fn::NAME: {                                                                                 \
    for (std::size_t l = 0; l < count; ++l) {                                                      \
      const T x = a0[l];                                                                           \
      const T y = a1[l];                                                                           \
      dst[l] = (EXPR);                                                                             \
    }                                                                                              \
    break;                                                                                         \
  }
#define SEISSOL_EXPR_PW_LOOP_TERNARY(NAME, EXPR)                                                   \
  case Fn::NAME: {                                                                                 \
    for (std::size_t l = 0; l < count; ++l) {                                                      \
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
      break;
    }
    }
  }

  // Stores run after every instruction, which is what makes state updates a
  // parallel assignment rather than a sequential one.
  for (const Store& store : stage.outputs) {
    const T* const __restrict src = slotPtr(store.source);
    T* const __restrict out = outputTile + static_cast<std::size_t>(store.target) * count;
    for (std::size_t l = 0; l < count; ++l) {
      out[l] = src[l];
    }
  }
  for (const Store& store : stage.persistent) {
    const T* const __restrict src = slotPtr(store.source);
    T* const __restrict out =
        persistent + static_cast<std::size_t>(store.target) * numPoints + first;
    for (std::size_t l = 0; l < count; ++l) {
      out[l] = src[l];
    }
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
      runStage(lowered_->precompute(),
               inputTile_.data(),
               outputTile_.data(),
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
      runStage(lowered_->run(),
               inputTile_.data(),
               outputTile_.data(),
               persistent,
               numPoints,
               first,
               count);
      io.scatter(first, count, outputTile_.data());
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

template class TileInterpreter<double>;
template class TileInterpreter<float>;
template void initialiseState<double>(const Program&, double*, std::size_t);
template void initialiseState<float>(const Program&, float*, std::size_t);

} // namespace seissol::expr
