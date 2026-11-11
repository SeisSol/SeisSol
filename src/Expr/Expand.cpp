// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/Expand.h"

#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Expr/Types.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/Interpolation.h"
#include "utils/logger.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace seissol::expr {

namespace {

namespace df = reader::datafield;

constexpr ValueType F64 = ValueType::F64;
constexpr ValueType I64 = ValueType::I64;
constexpr ValueType Bool = ValueType::Bool;

/// How far, in cells, the stencil of a time may reach beyond the resident window of a
/// time-dependent grid before the window counts as not serving it: rounding of the times, never a
/// slice it lacks.
constexpr double WindowTolerance = 1e-3;

[[noreturn]] void fail(const std::string& what) { throw std::invalid_argument("expr: " + what); }

/// One axis of a lookup, located: the first node of the stencil, the weights of its nodes, the last
/// node of the axis and the elements between two nodes (whole numbers); whether the coordinate lies
/// beyond the volume (Bool); and, on the time axis of a time-dependent grid, whether it lies beyond
/// the resident window where the file goes on (Bool).
struct Axis {
  NodeId start{NoNode};
  std::vector<NodeId> weights;
  NodeId top{NoNode};
  NodeId stride{NoNode};
  NodeId clamped{NoNode};
  NodeId stale{NoNode};
};

class Expander {
  public:
  explicit Expander(const Program& in) : in_(in), from_(in.arena()) {}

  ExpandedProgram expand() {
    Program& out = result_.program;
    out.setComputeType(in_.computeType());
    for (const auto& grid : in_.grids()) {
      static_cast<void>(out.internGrid(grid));
    }
    for (const auto& matrix : in_.matrices()) {
      static_cast<void>(out.internMatrix(matrix.name, matrix.shape));
    }
    for (const auto& block : in_.blocks()) {
      static_cast<void>(out.internBlock(block.name, block.length));
    }
    for (const auto& input : in_.inputs()) {
      out.addInput(input.name, input.type);
    }
    result_.uses.assign(in_.grids().size(), GridUse{});

    // what the outputs and the states need, made in the order of the arena, which is topological
    std::vector<char> reachable(from_.size(), 0);
    std::vector<NodeId> stack(in_.roots().begin(), in_.roots().end());
    for (const auto& state : in_.state()) {
      stack.push_back(state.root);
    }
    std::vector<NodeId> kids;
    while (!stack.empty()) {
      const NodeId id = stack.back();
      stack.pop_back();
      if (id == NoNode || reachable[static_cast<std::size_t>(id)] != 0) {
        continue;
      }
      reachable[static_cast<std::size_t>(id)] = 1;
      from_.children(id, kids);
      stack.insert(stack.end(), kids.begin(), kids.end());
    }
    mapped_.assign(from_.size(), NoNode);
    for (std::size_t index = 0; index < from_.size(); ++index) {
      if (reachable[index] != 0) {
        mapped_[index] = make(static_cast<NodeId>(index));
      }
    }

    for (std::size_t i = 0; i < in_.outputs().size(); ++i) {
      out.addOutput(in_.outputs()[i].name, in_.outputs()[i].type, mapped(in_.roots()[i]));
    }
    for (const auto& state : in_.state()) {
      out.addState(state.name, state.initial, mapped(state.root), state.type);
    }
    return std::move(result_);
  }

  private:
  [[nodiscard]] Arena& arena() { return result_.program.arena(); }
  [[nodiscard]] NodeId mapped(NodeId id) const { return mapped_[static_cast<std::size_t>(id)]; }

  // the numbers of the samples, and their arithmetic, each operation rounded on its own
  NodeId k(double value) { return arena().konst(value, F64); }
  NodeId sadd(NodeId x, NodeId y) { return arena().strictPw(Fn::Add, x, y, F64); }
  NodeId ssub(NodeId x, NodeId y) { return arena().strictPw(Fn::Sub, x, y, F64); }
  NodeId smul(NodeId x, NodeId y) { return arena().strictPw(Fn::Mul, x, y, F64); }
  // the whole numbers: indices, node counts, offsets
  NodeId n(std::int64_t value) { return arena().konst(static_cast<double>(value), I64); }
  NodeId iadd(NodeId x, NodeId y) { return arena().pw(Fn::Add, x, y, I64); }
  NodeId isub(NodeId x, NodeId y) { return arena().pw(Fn::Sub, x, y, I64); }
  NodeId imul(NodeId x, NodeId y) { return arena().pw(Fn::Mul, x, y, I64); }
  NodeId imin(NodeId x, NodeId y) { return arena().pw(Fn::Min, x, y, I64); }
  NodeId whole(NodeId x) { return arena().cast(I64, x); }
  NodeId real(NodeId x) { return arena().cast(F64, x); }
  // truth values: comparisons of numbers of either kind, and what combines them
  NodeId truth(bool value) { return arena().konst(value ? 1.0 : 0.0, Bool); }
  NodeId lt(NodeId x, NodeId y) { return arena().pw(Fn::Lt, x, y, Bool); }
  NodeId le(NodeId x, NodeId y) { return arena().pw(Fn::Le, x, y, Bool); }
  NodeId eq(NodeId x, NodeId y) { return arena().pw(Fn::Eq, x, y, Bool); }
  NodeId both(NodeId x, NodeId y) { return arena().pw(Fn::And, x, y, Bool); }
  NodeId either(NodeId x, NodeId y) { return arena().pw(Fn::Or, x, y, Bool); }
  /// `x` where `condition` (Bool), else `y`, in the type `type` of both
  NodeId choose(NodeId condition, NodeId x, NodeId y, ValueType type) {
    return arena().pw(Fn::Select, condition, x, y, type);
  }

  NodeId make(NodeId id) {
    const Node& node = from_[id];
    Arena& a = arena();
    switch (node.kind) {
    case Kind::Const:
      return a.konst(node.value, node.type);
    case Kind::Field:
      return a.field(a.channel(from_.channelName(node.ch)), node.type);
    case Kind::PW:
      if (node.strict) {
        return a.strictPw(node.fn, mapped(node.a), mapped(node.b), node.type);
      }
      switch (arity(node.fn)) {
      case 1:
        return a.pw(node.fn, mapped(node.a), node.type);
      case 2:
        return a.pw(node.fn, mapped(node.a), mapped(node.b), node.type);
      default:
        return a.pw(node.fn, mapped(node.a), mapped(node.b), mapped(node.c), node.type);
      }
    case Kind::Cast:
      return a.cast(node.type, mapped(node.a));
    case Kind::Contract:
      return a.contract(node.matrix, node.block, node.type);
    case Kind::Lookup:
      return lookup(node);
    case Kind::GridParam:
      return a.gridParam(node.grid, static_cast<GridField>(node.comp), node.axis);
    case Kind::GridLoad:
      return a.gridLoad(node.grid, node.comp, mapped(node.a));
    case Kind::Dx:
      return a.dx(node.axis, mapped(node.a));
    case Kind::Cumint:
      return a.cumint(mapped(node.a));
    case Kind::Fold:
      return a.fold(node.red, mapped(node.a));
    case Kind::Sample:
      return a.sample(mapped(node.a));
    }
    fail(std::string("a node of the kind ") + name(node.kind) + " cannot be written out");
  }

  NodeId lookup(const Node& node) {
    Arena& a = arena();
    const GridId grid = node.grid;
    if (grid < 0 || static_cast<std::size_t>(grid) >= in_.grids().size()) {
      fail("a lookup of grid " + std::to_string(grid) + ", which the program does not declare");
    }
    const df::GridDesc& desc = in_.grids()[static_cast<std::size_t>(grid)];
    const std::int32_t dimensions = node.argCount;
    GridUse& use = result_.uses[static_cast<std::size_t>(grid)];
    if (use.dimensions != 0 && use.dimensions != dimensions) {
      fail("the grid " + desc.path + " is sampled with " + std::to_string(use.dimensions) +
           " and with " + std::to_string(dimensions) + " coordinates");
    }
    use.dimensions = dimensions;
    use.components = std::max(use.components, node.comp + 1);

    // the coordinates in binary64, which typed() gave the type of the lookup
    const bool wide = resolve(node.type, in_.computeType()) == ComputeType::F64;
    const bool fails = desc.boundary == df::BoundaryMode::Fail;
    const NodeId* args = from_.args(node);
    std::vector<Axis> axes;
    for (std::int32_t d = 0; d < dimensions; ++d) {
      const NodeId coordinate = mapped(args[d]);
      axes.push_back(locate(grid, d, wide ? coordinate : real(coordinate), desc, fails));
    }

    // the offsets of the elements without the component, so that the lookups of several
    // components at one point share them
    const NodeId sample = reduce(grid, node.comp, axes, 0, n(0));

    if (fails) {
      NodeId outside = axes.front().clamped;
      for (std::size_t d = 1; d < axes.size(); ++d) {
        outside = either(outside, axes[d].clamped);
      }
      result_.checks.push_back(GridCheck{grid, CheckKind::Outside, outside});
    }
    for (const Axis& axis : axes) {
      if (axis.stale != NoNode) {
        result_.checks.push_back(GridCheck{grid, CheckKind::Window, axis.stale});
      }
    }
    return wide ? sample : a.cast(node.type, sample);
  }

  /// The stencil of the scheme on axis `d` at the coordinate `x`: locateAxis in Interpolation.cpp.
  Axis locate(GridId grid, std::int32_t d, NodeId x, const df::GridDesc& desc, bool fails) {
    Arena& a = arena();
    Axis axis;
    const NodeId num = a.gridParam(grid, GridField::Num, d);
    axis.top = isub(num, n(1));
    axis.stride = a.gridParam(grid, GridField::Stride, d);
    const NodeId top = real(axis.top);
    const NodeId raw = smul(ssub(x, a.gridParam(grid, GridField::Min, d)),
                            a.gridParam(grid, GridField::DeltaInv, d));
    // `raw > 0`, not `raw >= 0`: NaN clamps onto the lower end
    const NodeId above = lt(top, raw);
    const NodeId clampedRaw = choose(lt(k(0), raw), choose(above, top, raw, F64), k(0), F64);

    // beyond the volume: below or above it, or NaN. On the time axis of a grid with a window, the
    // volume is that of the file, whose slices beyond the window it knows; a time beyond the window
    // where the file goes on is no sample beyond the volume, but one the window should have served.
    const NodeId below = lt(raw, k(0));
    const NodeId undefined = eq(eq(raw, raw), truth(false));
    const bool windowed = desc.timeAxis.has_value() && *desc.timeAxis == d;
    if (windowed) {
      const NodeId noneBefore = eq(a.gridParam(grid, GridField::Before, d), n(0));
      const NodeId noneAfter = eq(a.gridParam(grid, GridField::After, d), n(0));
      if (fails) {
        axis.clamped = either(undefined, either(both(below, noneBefore), both(above, noneAfter)));
      }
      // The window serves a time where it holds the stencil that the file gives it
      // (placeTimeWindow() places it so up to the next update): the nodes [base + offset, base +
      // offset + width - 1] from the base cell, or the nearest node. In the cells of the window,
      // from the node the stencil reaches back from to the one it reaches on to, up to rounding.
      const df::StencilKernel kernel = df::kernelOf(desc.interpolation);
      const double reachBack = kernel.roundToNearest ? -0.5 : -static_cast<double>(kernel.offset);
      const double reachOn = kernel.roundToNearest ? -0.5
                                                   : static_cast<double>(kernel.width) - 2.0 +
                                                         static_cast<double>(kernel.offset);
      const NodeId early = lt(raw, k(reachBack - WindowTolerance));
      const NodeId late = lt(a.pw(Fn::Add, top, k(WindowTolerance - reachOn), F64), raw);
      axis.stale = either(both(early, eq(noneBefore, truth(false))),
                          both(late, eq(noneAfter, truth(false))));
    } else if (fails) {
      axis.clamped = either(either(below, undefined), above);
    }

    if (desc.interpolation == df::Interpolation::Nearest) {
      // half away from zero, as std::lround
      axis.start = whole(a.pw(Fn::Round, clampedRaw, F64));
      axis.weights = {k(1)};
      return axis;
    }

    // the base cell is the last one with a right neighbour: a coordinate on the last node is s = 1
    // in the cell before
    const NodeId twoNodes = le(n(2), num);
    const NodeId maxBase = choose(twoNodes, isub(num, n(2)), n(0), I64);
    const NodeId lowerReal = a.pw(Fn::Floor, clampedRaw, F64);
    const NodeId lower = whole(lowerReal);
    const NodeId past = lt(maxBase, lower);
    const NodeId index = choose(past, maxBase, lower, I64);
    const NodeId s =
        choose(past, choose(twoNodes, k(1), k(0), F64), ssub(clampedRaw, lowerReal), F64);
    const NodeId oneMinusS = ssub(k(1), s);

    if (desc.interpolation == df::Interpolation::Linear) {
      axis.start = index;
      axis.weights = {oneMinusS, s};
      return axis;
    }

    // cubic: four nodes from the one before the base cell, moved inside the axis at its ends,
    // where the weights fall back to linear ones between the base and its neighbour
    const df::StencilKernel kernel = df::kernelOf(desc.interpolation);
    const auto width = static_cast<std::int64_t>(kernel.width);
    const NodeId fourNodes = le(n(width), num);
    const NodeId maxStart = isub(num, n(width));
    const NodeId first = choose(fourNodes, iadd(index, n(kernel.offset)), n(0), I64);
    const NodeId start = choose(
        lt(first, n(0)),
        n(0),
        choose(lt(maxStart, first), choose(lt(n(0), maxStart), maxStart, n(0), I64), first, I64),
        I64);
    const NodeId slot = isub(index, start);
    const NodeId full = both(fourNodes, eq(slot, n(-kernel.offset)));

    // Keys' cubic convolution in the factored form of interpolationWeights()
    const NodeId s2 = smul(s, s);
    const NodeId s3 = smul(s2, s);
    const std::vector<NodeId> cubic = {
        smul(k(0.5), ssub(sadd(a.pw(Fn::Neg, s3, F64), smul(k(2), s2)), s)),
        smul(k(0.5), sadd(ssub(smul(k(3), s3), smul(k(5), s2)), k(2))),
        smul(k(0.5), sadd(sadd(smul(k(-3), s3), smul(k(4), s2)), s)),
        smul(k(0.5), ssub(s3, s2))};
    axis.start = start;
    for (std::size_t j = 0; j < cubic.size(); ++j) {
      // linearFallbackWeights(): 1 - s at the base, s at its neighbour
      const auto node = static_cast<std::int64_t>(j);
      const NodeId fallback =
          choose(eq(slot, n(node)), oneMinusS, choose(eq(slot, n(node - 1)), s, k(0), F64), F64);
      axis.weights.push_back(choose(full, cubic[j], fallback, F64));
    }
    return axis;
  }

  /// The tensor product over the axes from `d` on, of the elements of component `component` from
  /// `offset` on: the sum over the last axis innermost, each sum from zero in the order of its
  /// nodes, as the reduction of sampleBatch() forms them.
  NodeId reduce(GridId grid,
                std::int32_t component,
                const std::vector<Axis>& axes,
                std::size_t d,
                NodeId offset) {
    Arena& a = arena();
    const Axis& axis = axes[d];
    NodeId sum = k(0);
    for (std::size_t i = 0; i < axis.weights.size(); ++i) {
      // only the upper end needs a clamp: the start is not negative
      const NodeId shifted =
          i == 0 ? axis.start : iadd(axis.start, n(static_cast<std::int64_t>(i)));
      const NodeId along = imul(imin(shifted, axis.top), axis.stride);
      const NodeId at = d == 0 ? along : iadd(offset, along);
      const NodeId value = d + 1 == axes.size() ? a.gridLoad(grid, component, at)
                                                : reduce(grid, component, axes, d + 1, at);
      sum = sadd(sum, smul(axis.weights[i], value));
    }
    return sum;
  }

  const Program& in_;
  const Arena& from_;
  ExpandedProgram result_;
  std::vector<NodeId> mapped_;
};

} // namespace

ExpandedProgram expandLookups(const Program& program) { return Expander(program).expand(); }

std::vector<GridUse> gridUses(const Program& program) {
  const Arena& arena = program.arena();
  std::vector<GridUse> uses(program.grids().size());
  std::vector<char> seen(arena.size(), 0);
  std::vector<NodeId> stack(program.roots().begin(), program.roots().end());
  for (const auto& state : program.state()) {
    stack.push_back(state.root);
  }
  std::vector<NodeId> kids;
  while (!stack.empty()) {
    const NodeId id = stack.back();
    stack.pop_back();
    if (id == NoNode || seen[static_cast<std::size_t>(id)] != 0) {
      continue;
    }
    seen[static_cast<std::size_t>(id)] = 1;
    const Node& node = arena[id];
    if (node.kind == Kind::Lookup && node.grid >= 0 &&
        static_cast<std::size_t>(node.grid) < uses.size()) {
      GridUse& use = uses[static_cast<std::size_t>(node.grid)];
      if (use.dimensions != 0 && use.dimensions != node.argCount) {
        fail("the grid " + program.grids()[static_cast<std::size_t>(node.grid)].path +
             " is sampled with " + std::to_string(use.dimensions) + " and with " +
             std::to_string(node.argCount) + " coordinates");
      }
      use.dimensions = node.argCount;
      use.components = std::max(use.components, node.comp + 1);
    }
    arena.children(id, kids);
    stack.insert(stack.end(), kids.begin(), kids.end());
  }
  return uses;
}

double gridField(const df::ArrayView& view, GridField field, std::int32_t axis) {
  const auto d = static_cast<std::size_t>(axis);
  switch (field) {
  case GridField::Min:
    return view.min[d];
  case GridField::DeltaInv:
    return view.deltaInv[d];
  case GridField::Num:
    return static_cast<double>(view.num[d]);
  case GridField::Stride:
    return static_cast<double>(view.stride[d]);
  case GridField::ComponentStride:
    return static_cast<double>(view.componentStride);
  case GridField::Before:
    return static_cast<double>(view.before[d]);
  case GridField::After:
    return static_cast<double>(view.after[d]);
  }
  return 0.0;
}

void reportCheck(CheckKind kind, const std::string& path) {
  switch (kind) {
  case CheckKind::Outside:
    logError() << "expr: a point sampled the grid" << path
               << "beyond its volume, where its boundary mode is to fail.";
    return;
  case CheckKind::Window:
    logError() << "expr: a point sampled the time-dependent grid" << path
               << "at a time its window does not hold, where the file goes on. The window holds"
               << "the slices of the times from one synchronization point to the next: a script"
               << "that samples the grid at other times (t + 1, say) needs all slices of the file"
               << "resident (residentslices in &datafield); otherwise, SeisSol did not move the"
               << "window in time.";
    return;
  }
}

} // namespace seissol::expr
