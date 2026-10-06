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
#include "Expr/Program.h"

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

/// How a target spells the operands of a contraction, and what it has to know
/// about them: the form of every matrix, and the element types the binding
/// stores the matrices and blocks in -- those are baked into the source like the
/// column types of a device kernel.
struct ContractAddressing {
  std::vector<MatrixShape> shapes;      // by MatrixId
  std::vector<std::string> matrixTypes; // element type names, by MatrixId
  std::vector<std::string> blockTypes;  // element type names, by BlockId
  /// The point index of the lane at hand.
  std::string point;
  /// The dialect's unsigned 64-bit integer type.
  std::string indexType;
  /// The address-space qualifier of global memory ("__global " in OpenCL C, empty elsewhere).
  std::string global;
  /// Expressions for the operands: a pointer to the matrix, a pointer to the block data, the
  /// block's cell and coefficient strides in bytes, and its cell index (a pointer to unsigned
  /// int that may be null at run time).
  std::string (*matrix)(std::int32_t matrix){nullptr};
  std::string (*blockBase)(std::int32_t block){nullptr};
  std::string (*cellStride)(std::int32_t block){nullptr};
  std::string (*modeStride)(std::int32_t block){nullptr};
  std::string (*cellIndex)(std::int32_t block){nullptr};
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
  /// The COMPLETE store statement, without the trailing semicolon, given the
  /// name of the local holding the value.
  ///
  /// A statement rather than an assignment target, because not every dialect
  /// can produce one: OpenCL C has no operator overloading, so a device backend
  /// that wanted `store(...) = v` would need a second store form. Handing the
  /// value in lets each target spell the store however it can.
  std::string (*storeOutput)(std::int32_t index, const std::string& value){nullptr};
  std::string (*storePersistent)(std::int32_t slot, const std::string& value){nullptr};
  /// Null for a target that cannot read blocks; a contraction then throws.
  const ContractAddressing* contract{nullptr};
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
/// A contraction is emitted as the same loop the interpreter runs -- ascending
/// over the coefficients, one multiply-add after the other -- so a target that
/// compiles without contraction of floating-point expressions agrees with the
/// interpreter bit for bit.
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

/// True when either stage contracts a block.
[[nodiscard]] bool containsContraction(const LoweredProgram& lowered);

/// The C name of an element type, for the types the binding may store data in.
[[nodiscard]] const char* elementTypeName(reader::scripting::DataType type);

} // namespace seissol::expr::codegen

#endif // SEISSOL_SRC_EXPR_CODEGEN_H_
