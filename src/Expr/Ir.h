// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_IR_H_
#define SEISSOL_SRC_EXPR_IR_H_

// Shared expression IR. Derived from the derived-output `sderiv` IR (hash-consed
// Arena, Node, Fn/Red enums) with the additions the scripting path needs.
//
// DELTAS AGAINST THE EXISTING sderiv/ir.hpp — each one is a mechanical but
// non-local change to the existing port; see the porting notes at the bottom.
//   1. `arity(Fn)` is a table, not `f <= Fn::Sign`. The Fn enum is no longer
//      split in half, because Select is ternary and the comparison ops sit
//      between the old unary and binary blocks.
//   2. `Node` gains a third fixed child `c` (Select) and an (argBegin,argCount)
//      span for variable-arity nodes (Lookup coordinates, d = 1..6).
//   3. New `Kind::Lookup`: sample an external data grid. This is NOT the same
//      thing as `Kind::Sample`, which stays what it was in sderiv (materialize
//      at the output cadence). The names are close enough to be dangerous —
//      hence the split of the *word* "field": in this IR `Kind::Field` is a
//      named input channel, and an ASAGI-style data grid is a "grid" everywhere.
//   4. Const interning keys on the bit pattern, not on `==` (see NodeEq).

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

namespace seissol::expr {

using NodeId = std::int32_t;
inline constexpr NodeId NoNode = -1;

using GridId = std::int32_t;
inline constexpr GridId NoGrid = -1;

enum class Kind : std::uint8_t {
  Const,  // literal
  Field,  // named input channel (vx, x, t, group, …) — resolved via Binding
  PW,     // pointwise op, arity from the Fn
  Lookup, // sample an external data grid at the given coordinates  [NEW]
  Dx,     // physical-space derivative                       [derived-output]
  Cumint, // cumulative time integral                        [derived-output]
  Fold,   // temporal reducer                                [derived-output]
  Sample  // materialize at output cadence                   [derived-output]
};

// Pointwise ops. Arity is fixed per op; the enum order carries no meaning any
// more (see delta 1) — always go through arity().
enum class Fn : std::uint8_t {
  // unary
  Neg,
  Sqrt,
  Abs,
  Exp,
  Log,
  Log2,
  Log10,
  Sign,
  Floor,
  Ceil,
  Round,
  Rcp,
  Sin,
  Cos,
  Tan,
  Asin,
  Acos,
  Atan,
  Sinh,
  Cosh,
  Tanh,
  Asinh,
  Acosh,
  Atanh,
  Erf,
  // binary
  Add,
  Sub,
  Mul,
  Div,
  Pow,
  Min,
  Max,
  Atan2,
  Mod,
  // binary, boolean-valued (0.0 / 1.0 in the compute type)
  Lt,
  Le,
  Eq,
  And,
  Or,
  // ternary
  Select
};

constexpr int arity(Fn f) {
  switch (f) {
  case Fn::Neg:
  case Fn::Sqrt:
  case Fn::Abs:
  case Fn::Exp:
  case Fn::Log:
  case Fn::Log2:
  case Fn::Log10:
  case Fn::Sign:
  case Fn::Floor:
  case Fn::Ceil:
  case Fn::Round:
  case Fn::Rcp:
  case Fn::Sin:
  case Fn::Cos:
  case Fn::Tan:
  case Fn::Asin:
  case Fn::Acos:
  case Fn::Atan:
  case Fn::Sinh:
  case Fn::Cosh:
  case Fn::Tanh:
  case Fn::Asinh:
  case Fn::Acosh:
  case Fn::Atanh:
  case Fn::Erf:
    return 1;
  case Fn::Select:
    return 3;
  default:
    return 2;
  }
}

// Temporal reducers — unchanged from sderiv, unused by the scripting path.
enum class Red : std::uint8_t { Max, Min, Int, Sum, Mean, Rms, ArgMax, ArgMin, Last };

struct Node {
  Kind kind{};
  Fn fn{};                  // PW
  Red red{};                // Fold
  std::int32_t axis{};      // Dx
  std::int32_t ch{};        // Field: interned channel id
  GridId grid{NoGrid};      // Lookup
  std::int32_t comp{};      // Lookup: component within the grid
  double value{};           // Const
  NodeId a{NoNode};         // PW arg0 / Dx.x / Cumint.x / Fold.x / Sample.x
  NodeId b{NoNode};         // PW arg1
  NodeId c{NoNode};         // PW arg2 (Select only)
  std::int32_t argBegin{0}; // Lookup: coordinate children, span into Arena::args_
  std::int32_t argCount{0};
};

class Arena;

// Structural equality / hash, keyed per kind. Const compares BIT PATTERNS: the
// sderiv version used `x.value == y.value`, which makes the relation
// non-reflexive for NaN (a NaN key can never be found again in the pool) and
// leaves the 0.0/-0.0 hash-vs-equal agreement up to the standard library.
struct NodeEq {
  const Arena* arena{nullptr};
  bool operator()(const Node& x, const Node& y) const;
};
struct NodeHash {
  const Arena* arena{nullptr};
  std::size_t operator()(const Node& n) const;
};

// The interning arena. Owns every node and hands out stable NodeIds; equal
// subtrees get equal ids, so CSE falls out of construction and both the plan
// builder and the linearizer keep working on plain id equality.
class Arena {
  public:
  Arena();

  [[nodiscard]] const Node& operator[](NodeId id) const { return nodes_[id]; }
  [[nodiscard]] std::size_t size() const { return nodes_.size(); }

  [[nodiscard]] int channel(const std::string& name);
  [[nodiscard]] const std::string& channelName(int id) const { return channelNames_[id]; }
  [[nodiscard]] std::size_t channelCount() const { return channelNames_.size(); }

  NodeId konst(double v);
  NodeId field(int ch);
  NodeId field(const std::string& name) { return field(channel(name)); }
  NodeId pw(Fn f, NodeId x);
  NodeId pw(Fn f, NodeId x, NodeId y);
  NodeId pw(Fn f, NodeId x, NodeId y, NodeId z);
  NodeId lookup(GridId grid, std::int32_t component, const std::vector<NodeId>& coords);
  NodeId dx(int axis, NodeId x);
  NodeId cumint(NodeId x);
  NodeId fold(Red r, NodeId x);
  NodeId sample(NodeId x);

  // Ordered children, uniform across kinds. Returns {} for leaves.
  [[nodiscard]] std::vector<NodeId> children(NodeId id) const;
  // Allocation-free variant for hot walks; `out` is cleared first.
  void children(NodeId id, std::vector<NodeId>& out) const;

  [[nodiscard]] const NodeId* args(const Node& n) const { return args_.data() + n.argBegin; }

  private:
  NodeId intern(const Node& n);

  std::vector<Node> nodes_;
  std::vector<NodeId> args_;
  std::unordered_map<Node, NodeId, NodeHash, NodeEq> pool_;
  std::vector<std::string> channelNames_;
  std::unordered_map<std::string, int> channelIds_;
};

const char* name(Fn f);
const char* name(Red r);
bool fnFromName(const std::string& s, Fn& out);
bool redFromName(const std::string& s, Red& out);

} // namespace seissol::expr

// ---------------------------------------------------------------------------
// PORTING NOTES for the existing sderiv sources (all mechanical):
//
//   ir.hpp        replaced by this header. `is_unary(f)` disappears; call sites
//                 become `arity(f) == 1`.
//   lower.cpp     no change beyond the namespace and `is_unary`. `Kind::Lookup`
//                 must be added to the switches; treat it exactly like a leaf
//                 for cadence purposes if a grid is time-independent, and like
//                 a `t`-dependent node otherwise.
//   linearize.cpp add `lookup` and `select` instructions; the existing `memo`
//                 already handles sharing correctly.
//   codegen.hpp   Codegen::cx() is RECURSIVE WITHOUT MEMOISATION. On a DAG with
//                 shared subexpressions it re-expands each shared node inline at
//                 every reference, so the emitted source is exponential in the
//                 depth of the sharing (a diamond of depth d costs 2^d). This is
//                 already latent in the current code — linearize.cpp gets it
//                 right via `memo`, codegen does not. It must be fixed before
//                 the emitter is shared: emit SSA temporaries in linearized
//                 order instead of one nested expression.
//   interp.cpp    add the new Fn cases; the LANES/omp-simd structure is unchanged.
// ---------------------------------------------------------------------------

#endif // SEISSOL_SRC_EXPR_IR_H_
