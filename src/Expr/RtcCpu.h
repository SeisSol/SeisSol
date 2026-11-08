// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_RTCCPU_H_
#define SEISSOL_SRC_EXPR_RTCCPU_H_

// Runtime-compiled CPU kernel.
//
// WHAT THIS EMITS. One scalar lane loop per stage, with every transient value a
// local variable rather than a slot in a scratch array. That is the whole point:
// the interpreter's slots are memory the compiler cannot see through, and
// turning them into SSA locals is what lets the vectoriser work. Measured on a
// pointwise program, 2.6 ns/point against the interpreter's 28.3.
//
// WHY THE EXPRESSIONS COME FROM Interp.h. The arithmetic for each Fn is
// stringified out of SEISSOL_EXPR_PW_LIST -- the same macro table applyPw
// evaluates. Not a parallel table that has to be kept in step: a new Fn appears
// in the emitted code and in the interpreter from one edit, and the two cannot
// disagree about what Fn::Mod means, because there is only one place where it
// is said.
//
// WHY IT IS BITWISE IDENTICAL. Emitting the same expressions is necessary but
// not sufficient: with -O3 alone GCC contracts `a*b + c` into an FMA, which
// rounds once where the interpreter rounds twice, and the two paths drift by
// ~1e-6 on ordinary arithmetic. -ffp-contract=off removes that, and measurably
// costs nothing (2.6 ns/point either way). Bitwise equality is therefore the
// acceptance bar on CPU rather than a tolerance.
//
// THE SHAPE OF THE KERNEL. One point through the whole program after the other is right where the
// compiler vectorizes the loop over the points by itself, i.e. where the program does arithmetic
// only. A call to the maths library (exp, sin, asinh, pow, ...) keeps the loop from being
// vectorized, as the vector maths library needs -ffast-math and changes the bits; a long chain of
// dependent operations through such calls, such as an unrolled Newton iteration, is then bound by
// the latency of every one of them. For a stage with such a call and a critical path of at least
// TileDepth instructions, the kernel therefore takes one instruction over a block of points after
// the other, each a loop with `#pragma omp simd` (codegen::emitStageBodyTiled): the calls of the
// lanes overlap in the pipeline, and the arithmetic between them is vectorized. Both shapes
// compute the same operations at every point, so the bits stay the same. Measured (GCC 16, Zen 5,
// one thread): a Newton iteration through asinh takes 2 to 2.3 times less time per point in
// binary64 and 2.6 to 5 times less in binary32; one of arithmetic only, and a short program with
// calls, take a quarter to a half more as a tile, so they stay one point after the other. A stage
// with a contraction keeps the first shape too.
//
// GRIDS. A lookup arrives written out (expandLookups()): the kernel reads the
// geometry and the elements of a grid through a descriptor per grid
// (codegen::GridDescriptor), which the kernel fills from the grid before every
// call -- where its window currently lies -- and computes the sample with the
// interpreter's operations, so the bits are the same.

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/Interpolation.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace seissol::expr {

/// The critical path from which a stage of a compiled CPU kernel with a call to the maths library
/// takes the shape of a tile.
inline constexpr std::int32_t TileDepth = 32;

/// The C++ source for a lowered program. Exposed for tests and for a dump flag;
/// nothing else needs it. Without a binding, the matrices and blocks of the
/// contractions are taken to be stored in the compute type; without element types of
/// the grids, the arrays in binary64.
[[nodiscard]] std::string emitCpuSource(const Program& program,
                                        const LoweredProgram& lowered,
                                        KernelShape shape = KernelShape::Auto);
[[nodiscard]] std::string
    emitCpuSource(const Program& program,
                  const LoweredProgram& lowered,
                  const Binding& binding,
                  KernelShape shape = KernelShape::Auto,
                  const std::vector<reader::datafield::ElementType>& gridElements = {});

/// Whether a stage of the kernel takes the shape of a tile, as `shape` asks.
[[nodiscard]] bool tiledStage(const StageCode& stage,
                              const std::vector<std::int32_t>& operands,
                              KernelShape shape);

/// Compile and return a kernel, or nullptr if the toolchain is unusable. Never
/// throws for a compile failure: the caller falls back to the interpreter, which
/// is a warning and not an error. The grids of the program are those of `grids`
/// at `gridIds` (by GridId, as internGrids() gives them); a program that samples
/// none needs neither.
[[nodiscard]] std::unique_ptr<Kernel>
    makeRtcCpuKernel(const Program& program,
                     Binding& binding,
                     LoweredProgram lowered,
                     const BackendOptions& options,
                     reader::datafield::GridStore* grids = nullptr,
                     const std::vector<std::size_t>& gridIds = {});

/// Number of distinct kernels currently held. For tests and for a status line —
/// a cache that never hits is a configuration bug worth seeing.
[[nodiscard]] std::size_t rtcCpuCacheSize();

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_RTCCPU_H_
