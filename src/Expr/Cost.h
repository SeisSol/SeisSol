// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_COST_H_
#define SEISSOL_SRC_EXPR_COST_H_

// What a program costs per point, counted rather than measured.
//
// This is exact, not an estimate, and that is the point: after lowering there is
// no branching left and no data-dependent control flow, so the instruction list
// IS the execution. Every point runs the same sequence. A Select evaluates both
// arms because Interp.h says it does, so even that costs a fixed amount.
//
// WHAT IT IS FOR. Three things that were all guesswork before:
//   * deciding whether compiling is worth its 40 ms -- a program of twelve adds
//     will not repay it, one with four transcendentals per point will;
//   * putting a number next to a backend choice in the log, so "the interpreter
//     was slow" becomes a statement someone can check;
//   * sizing a batch, on the GPU in particular, where the useful question is
//     arithmetic intensity rather than raw flops -- which needs the byte counts
//     next to the operation counts, so they are here too.
//
// WHAT IT IS NOT. Not a timing model. A transcendental is one operation here and
// twenty to fifty cycles in reality, and the ratio moves between libm versions,
// between f32 and f64, and between host and device. The weighted figure carries
// a stated multiplier so the assumption is visible; the raw counts are what to
// compare across programs.

#include "Expr/Lower.h"
#include "Expr/Program.h"

#include <cstdint>
#include <string>

namespace seissol::expr {

/// Operation counts for ONE point, split by what actually differs in cost.
struct StageCost {
  /// Add, subtract, negate, and the comparison ops -- everything that is a
  /// single cheap instruction.
  std::uint64_t additive{0};
  /// Multiply. Separate from additive because a fused multiply-add is one
  /// instruction and the split is what tells you whether that applies.
  std::uint64_t multiplicative{0};
  /// Divide and reciprocal. An order of magnitude above a multiply and not
  /// pipelined the same way.
  std::uint64_t divisions{0};
  /// sqrt, exp, log, the trigonometric and hyperbolic family, erf, pow.
  std::uint64_t transcendentals{0};
  /// Grid samples. Counted separately because their cost is a stencil gather
  /// and a tensor-product reduction, not an instruction -- see Interpolation.h.
  std::uint64_t lookups{0};
  /// Loads from the gathered tile or the persistent buffer, and stores back.
  std::uint64_t loads{0};
  std::uint64_t stores{0};

  [[nodiscard]] std::uint64_t operations() const {
    return additive + multiplicative + divisions + transcendentals;
  }

  /// Operations with the expensive ones weighted, for a single comparable
  /// number. The multipliers are stated rather than tuned: a divide counts as
  /// `DivisionWeight` and a transcendental as `TranscendentalWeight`, both
  /// order-of-magnitude figures for a modern core in f64.
  static constexpr std::uint64_t DivisionWeight = 8;
  static constexpr std::uint64_t TranscendentalWeight = 20;

  [[nodiscard]] std::uint64_t weighted() const {
    return additive + multiplicative + DivisionWeight * divisions +
           TranscendentalWeight * transcendentals;
  }

  /// Bytes moved per point through the tile and the persistent buffer, at the
  /// program's compute type. Grid samples are excluded: what a lookup touches
  /// depends on the stencil and on cache residency, and a number that pretends
  /// otherwise is worse than none.
  [[nodiscard]] std::uint64_t bytes(ComputeType type) const;

  /// Weighted operations per byte. The figure that decides whether a device
  /// kernel will be bound by arithmetic or by memory.
  [[nodiscard]] double intensity(ComputeType type) const;
};

struct ProgramCost {
  StageCost precompute;
  StageCost run;

  /// One line for a log: counts, the weighted figure and the intensity.
  [[nodiscard]] std::string summary(ComputeType type) const;
};

/// Count what one point costs. Exact for everything except a lookup, whose cost
/// is reported as a count rather than folded in.
[[nodiscard]] ProgramCost cost(const LoweredProgram& lowered, ComputeType type);

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_COST_H_
