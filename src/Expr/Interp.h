// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_INTERP_H_
#define SEISSOL_SRC_EXPR_INTERP_H_

// Tiled batch/SoA interpreter over a LoweredProgram. Always available, and the
// correctness oracle every other backend is checked against — which is why the
// pointwise semantics live here, once, rather than being restated per backend.
//
// Dispatch is per instruction per tile, not per point: one switch amortised over
// `tileSize` lanes, with the lane loop left in a shape the vectoriser can take.
// That is the same structure as the sderiv interp.cpp, minus its text program
// format — the instruction stream arrives already lowered, already slot-
// allocated, and already split into stages.
//
// SEAMS. Two things this pass needs do not exist yet, so both enter through a
// narrow interface rather than a #include:
//   * GridSampler is Kind::Lookup's landing site and collapses into
//     datafield::Grid::sampleBatch when Paket 3 lands.
//   * TileIo is exactly Binding::gather / Binding::scatter, including the
//     dst[channel * count + lane] layout. Paket 4 supplies a two-line adapter
//     holding a (Binding, DataTable) pair.

#include "Expr/Ir.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace seissol::expr {

// --- pointwise semantics -----------------------------------------------------
//
// Written once, as a list, and expanded twice: into the scalar applyPw() below
// and into the vectorised lane loops in Interp.cpp. Two hand-written copies of
// 40 operator definitions would drift, and the drift would show up as a backend
// that is subtly wrong rather than one that fails to build.
//
// Semantics worth stating rather than inferring:
//   * Select evaluates BOTH branches. It is an operator over already-computed
//     values, not a branch. That is why layered models are handled by
//     partitioning the point set (Binding::groupRanges) instead of by a select
//     chain — the chain would read every grid at every point.
//   * Min/Max are fmin/fmax, so a NaN operand loses against a number.
//   * Sign is (x>0)-(x<0): zero for both signed zeros and for NaN.
//   * Rcp is an exact 1/x, not a hardware reciprocal estimate. A backend that
//     substitutes an approximation is no longer bit-comparable against this
//     interpreter and has to say so.
//   * Eq is an exact floating-point comparison, so NaN != NaN. It exists for
//     integer-valued channels (group, fault tag), where it is exact and safe.
#define SEISSOL_EXPR_PW_LIST(UNARY, BINARY, TERNARY)                                               \
  UNARY(Neg, (-x))                                                                                 \
  UNARY(Sqrt, (std::sqrt(x)))                                                                      \
  UNARY(Abs, (std::abs(x)))                                                                        \
  UNARY(Exp, (std::exp(x)))                                                                        \
  UNARY(Log, (std::log(x)))                                                                        \
  UNARY(Log2, (std::log2(x)))                                                                      \
  UNARY(Log10, (std::log10(x)))                                                                    \
  UNARY(Sign, (static_cast<T>((x > T(0)) - (x < T(0)))))                                           \
  UNARY(Floor, (std::floor(x)))                                                                    \
  UNARY(Ceil, (std::ceil(x)))                                                                      \
  UNARY(Round, (std::round(x)))                                                                    \
  UNARY(Rcp, (T(1) / x))                                                                           \
  UNARY(Sin, (std::sin(x)))                                                                        \
  UNARY(Cos, (std::cos(x)))                                                                        \
  UNARY(Tan, (std::tan(x)))                                                                        \
  UNARY(Asin, (std::asin(x)))                                                                      \
  UNARY(Acos, (std::acos(x)))                                                                      \
  UNARY(Atan, (std::atan(x)))                                                                      \
  UNARY(Sinh, (std::sinh(x)))                                                                      \
  UNARY(Cosh, (std::cosh(x)))                                                                      \
  UNARY(Tanh, (std::tanh(x)))                                                                      \
  UNARY(Asinh, (std::asinh(x)))                                                                    \
  UNARY(Acosh, (std::acosh(x)))                                                                    \
  UNARY(Atanh, (std::atanh(x)))                                                                    \
  UNARY(Erf, (std::erf(x)))                                                                        \
  BINARY(Add, (x + y))                                                                             \
  BINARY(Sub, (x - y))                                                                             \
  BINARY(Mul, (x * y))                                                                             \
  BINARY(Div, (x / y))                                                                             \
  BINARY(Pow, (std::pow(x, y)))                                                                    \
  BINARY(Min, (std::fmin(x, y)))                                                                   \
  BINARY(Max, (std::fmax(x, y)))                                                                   \
  BINARY(Atan2, (std::atan2(x, y)))                                                                \
  BINARY(Mod, (std::fmod(x, y)))                                                                   \
  BINARY(Lt, (x < y ? T(1) : T(0)))                                                                \
  BINARY(Le, (x <= y ? T(1) : T(0)))                                                               \
  BINARY(Eq, (x == y ? T(1) : T(0)))                                                               \
  BINARY(And, ((x != T(0) && y != T(0)) ? T(1) : T(0)))                                            \
  BINARY(Or, ((x != T(0) || y != T(0)) ? T(1) : T(0)))                                             \
  TERNARY(Select, (x != T(0) ? y : z))

// Scalar evaluation of one pointwise op. `args` holds arity(f) values.
template <typename T>
T applyPw(Fn f, const T* args) {
  switch (f) {
#define SEISSOL_EXPR_PW_SCALAR_UNARY(NAME, EXPR)                                                   \
  case Fn::NAME: {                                                                                 \
    const T x = args[0];                                                                           \
    return (EXPR);                                                                                 \
  }
#define SEISSOL_EXPR_PW_SCALAR_BINARY(NAME, EXPR)                                                  \
  case Fn::NAME: {                                                                                 \
    const T x = args[0];                                                                           \
    const T y = args[1];                                                                           \
    return (EXPR);                                                                                 \
  }
#define SEISSOL_EXPR_PW_SCALAR_TERNARY(NAME, EXPR)                                                 \
  case Fn::NAME: {                                                                                 \
    const T x = args[0];                                                                           \
    const T y = args[1];                                                                           \
    const T z = args[2];                                                                           \
    return (EXPR);                                                                                 \
  }
    SEISSOL_EXPR_PW_LIST(
        SEISSOL_EXPR_PW_SCALAR_UNARY, SEISSOL_EXPR_PW_SCALAR_BINARY, SEISSOL_EXPR_PW_SCALAR_TERNARY)
#undef SEISSOL_EXPR_PW_SCALAR_UNARY
#undef SEISSOL_EXPR_PW_SCALAR_BINARY
#undef SEISSOL_EXPR_PW_SCALAR_TERNARY
  }
  return T(0);
}

// --- seams -------------------------------------------------------------------

// Kind::Lookup lowers onto this. Coordinates and results are SoA in both
// directions, matching the Grid.h contract: `coords[d]` is a contiguous run of
// `count` values for dimension d.
class GridSampler {
  public:
  GridSampler() = default;
  virtual ~GridSampler() = default;
  GridSampler(const GridSampler&) = delete;
  GridSampler& operator=(const GridSampler&) = delete;
  GridSampler(GridSampler&&) = delete;
  GridSampler& operator=(GridSampler&&) = delete;

  virtual void sampleBatch(GridId grid,
                           std::int32_t component,
                           const double* const* coords,
                           std::int32_t dimension,
                           std::size_t count,
                           double* dst) const = 0;
  virtual void sampleBatch(GridId grid,
                           std::int32_t component,
                           const float* const* coords,
                           std::int32_t dimension,
                           std::size_t count,
                           float* dst) const = 0;
};

// Binding::gather / Binding::scatter, with the table already bound in.
template <typename T>
class TileIo {
  public:
  TileIo() = default;
  virtual ~TileIo() = default;
  TileIo(const TileIo&) = delete;
  TileIo& operator=(const TileIo&) = delete;
  TileIo(TileIo&&) = delete;
  TileIo& operator=(TileIo&&) = delete;

  // dst[inputIndex * count + lane]
  virtual void gather(std::size_t first, std::size_t count, T* dst) const = 0;
  // src[outputIndex * count + lane]
  virtual void scatter(std::size_t first, std::size_t count, const T* src) = 0;
};

// A contiguous run of points the tiling must not cross. Paket 4 fills this from
// Binding::groupRanges(); a program without a group input gets one range.
struct PointRange {
  std::size_t begin{0};
  std::size_t end{0};
};

// --- tile sizing -------------------------------------------------------------

inline constexpr std::size_t DefaultTileBudgetBytes = 16384; // half a 32 KiB L1D
inline constexpr std::size_t TileLaneGranularity = 8;
inline constexpr std::size_t MinTileSize = 8;
inline constexpr std::size_t MaxTileSize = 2048;

// The tile width follows from the live-value count, not the other way round: the
// value tile is what has to stay in L1, and after slot reuse the live count is a
// property of the program. The budget covers the value tile only — the input and
// output tiles and any grid stencil working set share the same L1, which is why
// the default is half of it rather than all of it.
std::size_t chooseTileSize(std::int32_t peakSlots, ComputeType type, std::size_t budgetBytes);

struct InterpreterOptions {
  std::size_t tileSize{0}; // 0 -> chooseTileSize
  std::size_t budgetBytes{DefaultTileBudgetBytes};
};

// --- interpreter -------------------------------------------------------------

template <typename T>
class TileInterpreter {
  public:
  // `program` supplies the input/output counts and the compute type; `lowered`
  // the instruction stream. Both must outlive the interpreter.
  TileInterpreter(const Program& program,
                  const LoweredProgram& lowered,
                  const GridSampler* sampler,
                  const InterpreterOptions& options);

  [[nodiscard]] std::size_t tileSize() const { return tileSize_; }

  // Runs the Precompute stage over every point and fills the hoisted persistent
  // slots. Must happen once, before the first run(), and again after anything
  // that invalidates the invariants it was told about. Explicit rather than lazy
  // on purpose: the cost belongs in the profile where it is spent.
  void precompute(const TileIo<T>& io,
                  std::size_t numPoints,
                  T* persistent,
                  const std::vector<PointRange>& partitions = {});

  // Runs the Run stage: gather, evaluate, scatter, and write state back. State
  // updates are parallel assignment — every state root sees the values from
  // before this call, because all stores happen after all instructions.
  void run(TileIo<T>& io,
           std::size_t numPoints,
           T* persistent,
           const std::vector<PointRange>& partitions = {});

  private:
  void runStage(const StageCode& stage,
                const T* inputTile,
                T* outputTile,
                T* persistent,
                std::size_t numPoints,
                std::size_t first,
                std::size_t count);

  const LoweredProgram* lowered_;
  const GridSampler* sampler_;
  std::size_t tileSize_{0};
  std::vector<T> scratch_;
  std::vector<T> inputTile_;
  std::vector<T> outputTile_;
  std::int32_t numInputs_{0};
  std::int32_t numOutputs_{0};
};

extern template class TileInterpreter<double>;
extern template class TileInterpreter<float>;

// Fills the declared state slots from StateSpec::initial. Binding owns the
// buffer and calls this on bind; exposed here so the storage layout has exactly
// one definition. Persistent slots [0, state().size()) are the declared states.
template <typename T>
void initialiseState(const Program& program, T* persistent, std::size_t numPoints);

extern template void initialiseState<double>(const Program&, double*, std::size_t);
extern template void initialiseState<float>(const Program&, float*, std::size_t);

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_INTERP_H_
