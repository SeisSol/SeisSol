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
  /// has no <cmath> to pull `std` in from. For OpenCL C, whose built-ins are
  /// overloaded for float and double alike.
  Unqualified,
  /// The C names, which name their type: `sqrtf(...)` in binary32, `sqrt(...)` in binary64,
  /// `fabsf`/`fabs` for `abs`. For CUDA and HIP, whose unqualified names reach the float overloads
  /// only through the C++ headers of the toolchain at hand; where an older one lacks them, a call
  /// in binary32 would silently compute in binary64.
  CNames
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
  /// The persistent slots below this are the declared states, which loadState and storeState
  /// read and write; the others hold hoisted values, and loadPersistent and storePersistent get
  /// their index among those (the slot less stateSlots).
  std::int32_t stateSlots{0};
  /// Expression for `LoadPersistent` of a declared state.
  std::string (*loadState)(std::int32_t state){nullptr};
  /// Expression for `LoadPersistent` of hoisted value `i`.
  std::string (*loadPersistent)(std::int32_t slot){nullptr};
  /// The COMPLETE store statement, without the trailing semicolon, given the
  /// name of the local holding the value.
  ///
  /// A statement rather than an assignment target, because not every dialect
  /// can produce one: OpenCL C has no operator overloading, so a device backend
  /// that wanted `store(...) = v` would need a second store form. Handing the
  /// value in lets each target spell the store however it can.
  std::string (*storeOutput)(std::int32_t index, const std::string& value){nullptr};
  std::string (*storeState)(std::int32_t state, const std::string& value){nullptr};
  std::string (*storePersistent)(std::int32_t slot, const std::string& value){nullptr};
  /// Null for a target that cannot read blocks; a contraction then throws.
  const ContractAddressing* contract{nullptr};
  /// The parts of a written-out lookup (expandLookups()): a field of the geometry of grid `grid`
  /// as the descriptor holds it (Opcode::GridParam); the element of component `component` of its
  /// array at `offset`, an expression of type I64, in binary64 (Opcode::GridLoad); and the
  /// statement that reports a check of the kind `kind` where `value` is true (StageCode::checks).
  /// Null for a target that samples no grids; a lookup then throws.
  std::string (*gridParam)(GridId grid, GridField field, std::int32_t axis){nullptr};
  std::string (*gridLoad)(GridId grid, std::int32_t component, const std::string& offset){nullptr};
  std::string (*storeCheck)(GridId grid, CheckKind kind, const std::string& value){nullptr};
};

/// The descriptor of the array of a grid, through which a compiled kernel reads it: the byte image
/// of `SeissolExprGrid` in the source of the kernel (emitGridDescriptor()). The CPU kernel reads
/// one in host memory, a device kernel one on the device, at an address that stays the same, so
/// that one recorded in a graph reads the window that is resident when it runs. Every field is
/// eight bytes wide, so the image has no padding to agree on.
struct GridDescriptor {
  /// The address of the first element, where the kernel runs.
  std::uint64_t data{0};
  double min[reader::datafield::MaxGridDimensions]{};
  double deltaInv[reader::datafield::MaxGridDimensions]{};
  std::uint64_t num[reader::datafield::MaxGridDimensions]{};
  std::uint64_t stride[reader::datafield::MaxGridDimensions]{};
  std::uint64_t before[reader::datafield::MaxGridDimensions]{};
  std::uint64_t after[reader::datafield::MaxGridDimensions]{};
  std::uint64_t componentStride{1};
  /// Of an axis of knots, the addresses of its knots and of the reciprocals of their spacings,
  /// where the kernel runs (ArrayView::knots, ArrayView::inverseSpacings); zero for a uniform axis.
  std::uint64_t knots[reader::datafield::MaxGridDimensions]{};
  std::uint64_t inverseSpacings[reader::datafield::MaxGridDimensions]{};
  /// Set by a point that sampled a grid that fails beyond its volume beyond it
  /// (CheckKind::Outside), and by one whose time the window of a time-dependent grid did not serve
  /// (CheckKind::Window), which a kernel cannot report itself; whoever keeps the descriptor reads
  /// them and does.
  std::uint64_t outside{0};
  std::uint64_t stale{0};
};

/// The descriptor of the array `view`, with its elements at `data` -- where the kernel runs: on the
/// host, view.data.
[[nodiscard]] GridDescriptor gridDescriptor(const reader::datafield::ArrayView& view,
                                            const void* data);

/// The bytes of the array of `view` a kernel may read: up to its last element.
[[nodiscard]] std::size_t gridBytes(const reader::datafield::ArrayView& view);

/// The definition of `SeissolExprGrid` (GridDescriptor), with `u64` the dialect's unsigned 64-bit
/// integer type; and per grid `g` that `uses` samples (by GridId) the function
/// `seissol_expr_load<g>(grid, offset)`, which reads the element at `offset` of its array, of the
/// type in `elements`, as binary64. `qualifier` comes before a function (`__device__ inline`, say),
/// `global` before a pointer into global memory (`__global ` in OpenCL C).
void emitGridDescriptor(std::ostringstream& out,
                        const std::vector<reader::datafield::ElementType>& elements,
                        const std::vector<GridUse>& uses,
                        const std::string& u64,
                        const std::string& qualifier,
                        const std::string& global);

/// The macros of the arithmetic that is computed as written (Instruction::strict), with the
/// definitions of a target, unless defined before: SEISSOL_EXPR_DADD, SEISSOL_EXPR_DSUB and
/// SEISSOL_EXPR_DMUL, each of two operands in binary64.
void emitStrictMacros(std::ostringstream& out,
                      const std::string& add,
                      const std::string& sub,
                      const std::string& mul);

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
/// A lookup samples through the addressing (StageAddressing::lookup), at its coordinates converted
/// to binary64, as the interpreter does; it throws std::invalid_argument for a target without one
/// -- written-out lookups reach a target that does not sample grids only by a programming error
/// rather than an unsupported model.
///
/// Each instruction computes in its type (Instruction::type), each transient is a local of the
/// type of its slot; the addressing reads and writes the data in `dataType`
/// (LoweredProgram::dataType()), and the body converts between the two where they differ.
void emitStageBody(std::ostringstream& out,
                   const StageCode& stage,
                   const std::vector<std::int32_t>& operands,
                   ComputeType dataType,
                   MathStyle style,
                   const StageAddressing& addressing,
                   const char* indent);

/// The same body in the shape of a tile: the lanes of a block of `block` points go through one
/// instruction after the other, each a loop over the lanes with `#pragma omp simd`, the transients
/// arrays of the block. The caller supplies the loop over the points `l < count` no longer: the
/// body loops over the blocks of [0, count) itself, with `l` the lane in [0, count) where the
/// addressing reads it. The lanes of a long chain of dependent operations then overlap in the
/// vector units and in the pipeline, where one point after the other waits for each result; every
/// lane computes the same operations as before, so the bits stay the same. Throws on a contraction,
/// which reads a block per lane.
void emitStageBodyTiled(std::ostringstream& out,
                        const StageCode& stage,
                        const std::vector<std::int32_t>& operands,
                        ComputeType dataType,
                        MathStyle style,
                        const StageAddressing& addressing,
                        std::size_t block,
                        const char* indent);

/// Whether `fn` is a call to the maths library in the emitted code, which a compiler does not
/// vectorize without the vector maths library (and -ffast-math): exp, log, the trigonometric and
/// hyperbolic functions, erf, pow, atan2 and fmod.
[[nodiscard]] bool callsLibrary(Fn fn);

/// The longest chain of dependent instructions of `stage`.
[[nodiscard]] std::int32_t criticalPath(const StageCode& stage,
                                        const std::vector<std::int32_t>& operands);

/// True when either stage contracts a block.
[[nodiscard]] bool containsContraction(const LoweredProgram& lowered);

/// The C name of an element type, for the types the binding may store data in.
[[nodiscard]] const char* elementTypeName(reader::scripting::DataType type);

/// The C name of a compute type: float, double, SeissolI64 (emitIntegerType()) or bool.
[[nodiscard]] const char* typeName(ComputeType type);

/// The definition of SeissolI64, the whole numbers of the lowering (ComputeType::I64): `long long`,
/// or `long` in OpenCL C.
void emitIntegerType(std::ostringstream& out, bool openCl);

/// Writes `text` to the file `name` in the directory SEISSOL_EXPR_DUMP names, if it names one: the
/// sources the compiled backends generate, and what the device compilers make of them, for a
/// developer to read. Does nothing otherwise.
void dump(const std::string& name, const std::string& text);

} // namespace seissol::expr::codegen

#endif // SEISSOL_SRC_EXPR_CODEGEN_H_
