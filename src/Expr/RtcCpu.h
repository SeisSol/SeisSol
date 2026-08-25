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
// WHAT IT REFUSES. Programs containing Opcode::Lookup. A lookup is a batch call
// into GridSampler, not a lane-local expression, and threading the sampler
// through a compiled kernel is separate work. makeKernel falls back to the
// interpreter, loudly, so a model with a data field still runs.

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"

#include <cstdint>
#include <memory>
#include <string>

namespace seissol::expr {

/// The C++ source for a lowered program. Exposed for tests and for a dump flag;
/// nothing else needs it.
[[nodiscard]] std::string emitCpuSource(const Program& program, const LoweredProgram& lowered);

/// True when this lowering can be compiled at all — today, when it contains no
/// Opcode::Lookup.
[[nodiscard]] bool cpuCompilable(const LoweredProgram& lowered);

/// Compile and return a kernel, or nullptr if the program is not compilable or
/// the toolchain is unusable. Never throws for a compile failure: the caller
/// falls back to the interpreter, which is a warning and not an error.
[[nodiscard]] std::unique_ptr<Kernel> makeRtcCpuKernel(const Program& program,
                                                       Binding& binding,
                                                       LoweredProgram lowered,
                                                       const BackendOptions& options);

/// Number of distinct kernels currently held. For tests and for a status line —
/// a cache that never hits is a configuration bug worth seeing.
[[nodiscard]] std::size_t rtcCpuCacheSize();

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_RTCCPU_H_
