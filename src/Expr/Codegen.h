// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_CODEGEN_H_
#define SEISSOL_SRC_EXPR_CODEGEN_H_

// The parts of code generation that every compiled backend shares.
//
// There is exactly one place in this codebase that says what Fn::Mod computes,
// and it is SEISSOL_EXPR_PW_LIST in Interp.h. The interpreter evaluates that
// table; everything here stringifies it. A backend that spelled the arithmetic
// again would be a second definition that agrees with the first only by
// inspection -- and it is precisely the rarely-taken branch, the one nobody
// re-reads, where the two would drift.

#include "Expr/Ir.h"
#include "Expr/Lower.h"

#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

namespace seissol::expr::codegen {

/// How a target spells the standard maths functions.
enum class MathStyle : std::uint8_t {
  /// `std::sqrt(...)`. Host C++ with <cmath>.
  Namespaced,
  /// `sqrt(...)`. Device code, where the built-ins are unqualified and NVRTC
  /// has no <cmath> to pull `std` in from.
  Unqualified
};

/// Where an instruction reads its inputs from and writes its outputs to. The
/// arithmetic is identical across backends; only these differ, because a CPU
/// kernel works on a gathered tile and a device kernel reads the strided views
/// straight from global memory.
struct StageAddressing {
  /// Expression for `LoadInput` with slot `i`, e.g. "inputTile[3ul * count + l]".
  std::string (*loadInput)(std::int32_t index){nullptr};
  /// Expression for `LoadPersistent` with slot `i`.
  std::string (*loadPersistent)(std::int32_t slot){nullptr};
  /// Assignment target for an output store.
  std::string (*storeOutput)(std::int32_t index){nullptr};
  /// Assignment target for a persistent store.
  std::string (*storePersistent)(std::int32_t slot){nullptr};
};

/// The expression text for `fn`, straight out of the interpreter's table, with
/// `x`, `y`, `z` still as placeholders.
[[nodiscard]] const char* expressionText(Fn fn);

/// Substitute `x`, `y`, `z` and `T` as whole IDENTIFIERS, and optionally strip
/// the `std::` qualification.
///
/// Not a plain string replace, and the reason is worth stating: `std::exp(x)`
/// contains an `x` inside `exp`. A naive replace turns that into
/// `std::e<operand>p(<operand>)`, which still compiles for some operand names
/// and computes nonsense.
[[nodiscard]] std::string substitute(const std::string& text,
                                     const std::string& x,
                                     const std::string& y,
                                     const std::string& z,
                                     const std::string& computeType,
                                     MathStyle style);

/// Name of the local holding transient slot `slot`.
[[nodiscard]] std::string slotName(std::int32_t slot);

/// A literal with enough digits to round-trip an IEEE double, so a Const in the
/// emitted source is the bit pattern the interpreter materialises. Not
/// std::to_string, which truncates at six decimals and would quietly emit a
/// different program than the one that was lowered.
[[nodiscard]] std::string literal(double value, const std::string& computeType);

/// Emit one stage's body: the local declarations, the instructions in order,
/// and the stores. The caller supplies the loop header and footer, because that
/// is where a lane loop and a grid-stride loop differ.
///
/// Every transient becomes a local variable rather than a slot in a scratch
/// array. That is not cosmetic: the interpreter's slots are memory the compiler
/// cannot see through, and turning them into SSA locals is what lets the
/// vectoriser work at all.
///
/// Throws std::invalid_argument on Opcode::Lookup -- callers gate on
/// containsLookup() first, and reaching here with one is a programming error
/// rather than an unsupported model.
void emitStageBody(std::ostringstream& out,
                   const StageCode& stage,
                   const std::vector<std::int32_t>& operands,
                   const std::string& computeType,
                   MathStyle style,
                   const StageAddressing& addressing,
                   const char* indent);

/// True when either stage samples a grid. No compiled backend supports that
/// yet: a lookup is a batch call into GridSampler, not a lane-local expression.
[[nodiscard]] bool containsLookup(const LoweredProgram& lowered);

} // namespace seissol::expr::codegen

#endif // SEISSOL_SRC_EXPR_CODEGEN_H_
