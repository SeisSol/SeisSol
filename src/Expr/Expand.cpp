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
    for (const auto& check : in_.checks()) {
      stack.push_back(check.root);
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
    // a point fails a check where its condition is zero; a NaN is true, as every value but zero
    for (std::size_t i = 0; i < in_.checks().size(); ++i) {
      const NodeId condition = mapped(in_.checks()[i].root);
      const ValueType type = arena()[condition].type;
      out.addCheck(condition, in_.checks()[i].message);
      result_.checks.push_back(GridCheck{NoGrid,
                                         CheckKind::User,
                                         eq(condition, arena().konst(0.0, type)),
                                         static_cast<std::int32_t>(i)});
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
  NodeId imax(NodeId x, NodeId y) { return arena().pw(Fn::Max, x, y, I64); }
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
      fail("the grid " + desc.label() + " is sampled with " + std::to_string(use.dimensions) +
           " and with " + std::to_string(dimensions) + " coordinates");
    }
    use.dimensions = dimensions;
    use.components = std::max(use.components, node.comp + 1);

    // the coordinates in binary64, which typed() gave the type of the lookup
    const bool wide = resolve(node.type, in_.computeType()) == ComputeType::F64;
    const NodeId* args = from_.args(node);
    std::vector<Axis> axes;
    for (std::int32_t d = 0; d < dimensions; ++d) {
      const NodeId coordinate = mapped(args[d]);
      axes.push_back(locate(grid, d, wide ? coordinate : real(coordinate), desc));
    }

    // the checks before the sample, so that their parts end before it: beyond the axes that fail
    NodeId outside = NoNode;
    for (const Axis& axis : axes) {
      if (axis.clamped != NoNode) {
        outside = outside == NoNode ? axis.clamped : either(outside, axis.clamped);
      }
    }
    if (outside != NoNode) {
      addCheck(grid, CheckKind::Outside, outside);
    }
    for (const Axis& axis : axes) {
      if (axis.stale != NoNode) {
        addCheck(grid, CheckKind::Window, axis.stale);
      }
    }

    // the offsets of the elements without the component, so that the lookups of several
    // components at one point share them
    const NodeId sample = reduce(grid, node.comp, axes, 0, n(0));
    return wide ? sample : a.cast(node.type, sample);
  }

  /// Adds `root` to the check of `grid` and `kind`: one truth value per grid and kind for the whole
  /// program, the `or` of those of its lookups, made as each lookup is made -- so that a kernel
  /// keeps one value per check alive to its end, not one per lookup.
  void addCheck(GridId grid, CheckKind kind, NodeId root) {
    for (GridCheck& check : result_.checks) {
      if (check.grid == grid && check.kind == kind) {
        check.root = either(check.root, root);
        return;
      }
    }
    result_.checks.push_back(GridCheck{grid, kind, root});
  }

  /// The stencil of the scheme of axis `d` at the coordinate `x`: locateAxis in Interpolation.cpp
  /// for a uniform axis, locateKnots() for one of knots.
  Axis locate(GridId grid, std::int32_t d, NodeId x, const df::GridDesc& desc) {
    const auto axisIndex = static_cast<std::size_t>(d);
    const df::Interpolation interpolation = desc.interpolationOf(axisIndex);
    const bool fails = desc.boundaryOf(axisIndex) == df::BoundaryMode::Fail;
    if (desc.axisOf(axisIndex) == df::AxisKind::Knots) {
      return locateKnots(grid, d, x, desc, interpolation, fails);
    }
    Arena& a = arena();
    Axis axis;
    const NodeId num = a.gridParam(grid, GridField::Num, d);
    axis.top = isub(num, n(1));
    axis.stride = a.gridParam(grid, GridField::Stride, d);
    const NodeId top = real(axis.top);
    const NodeId raw = smul(ssub(x, a.gridParam(grid, GridField::Min, d)),
                            a.gridParam(grid, GridField::DeltaInv, d));
    // `raw > 0`, not `raw >= 0`: NaN clamps onto the lower end
    const NodeId positive = lt(k(0), raw);
    const NodeId above = lt(top, raw);
    const NodeId clampedRaw = choose(positive, choose(above, top, raw, F64), k(0), F64);

    // beyond the volume: below or above it, or NaN. On the time axis of a grid with a window, the
    // volume is that of the file, whose slices beyond the window it knows; a time beyond the window
    // where the file goes on is no sample beyond the volume, but one the window should have served.
    const NodeId below = lt(raw, k(0));
    const NodeId undefined = eq(eq(raw, raw), truth(false));
    const bool windowed = desc.timeAxis.has_value() && *desc.timeAxis == d;
    NodeId someBefore = NoNode;
    NodeId someAfter = NoNode;
    if (windowed) {
      someBefore = lt(n(0), a.gridParam(grid, GridField::Before, d));
      someAfter = lt(n(0), a.gridParam(grid, GridField::After, d));
      if (fails) {
        axis.clamped = either(undefined,
                              either(both(below, eq(someBefore, truth(false))),
                                     both(above, eq(someAfter, truth(false)))));
      }
    } else if (fails) {
      axis.clamped = either(either(below, undefined), above);
    }
    // Whether the window does not hold the stencil the file gives the time, where the file goes on:
    // `early` where the stencil reaches before the window, `late` where it reaches beyond it. The
    // window keeps a slice to spare on either side of the times it is placed for (TimeWindowSlack),
    // so the predicates need not be exact at its ends -- they also hold in its first cell, or on
    // the last node but one for the cubic stencil, which those times never reach -- and are those
    // the stencil computes anyway, in whole numbers where they are new: no binary64 operation (a
    // negated comparison would be one again on a device), and no tolerance.
    const auto stale = [&](NodeId early, NodeId late) {
      if (windowed) {
        axis.stale = either(both(early, someBefore), both(late, someAfter));
      }
    };

    if (interpolation == df::Interpolation::Nearest) {
      // half away from zero, as std::lround
      axis.start = whole(a.pw(Fn::Round, clampedRaw, F64));
      axis.weights = {k(1)};
      stale(eq(axis.start, n(0)), above);
      return axis;
    }
    if (interpolation == df::Interpolation::Previous || interpolation == df::Interpolation::Next) {
      // the nodes below the coordinate: those of an index below raw, all beyond the axis
      const NodeId below = choose(above, num, whole(a.pw(Fn::Ceil, clampedRaw, F64)), I64);
      axis.start = interpolation == df::Interpolation::Previous ? imax(isub(below, n(1)), n(0))
                                                                : imin(below, axis.top);
      axis.weights = {k(1)};
      stale(eq(axis.start, n(0)), above);
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

    if (interpolation == df::Interpolation::Linear) {
      axis.start = index;
      axis.weights = {oneMinusS, s};
      stale(eq(lower, n(0)), above);
      return axis;
    }

    // cubic: four nodes from the one before the base cell, moved inside the axis at its ends,
    // where the weights fall back to linear ones between the base and its neighbour
    const df::StencilKernel kernel = df::kernelOf(interpolation);
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
    // the stencil moved inside the window, at its start or at its end, or a window too short for it
    const NodeId cramped = eq(fourNodes, truth(false));
    stale(either(cramped, lt(slot, n(-kernel.offset))),
          either(cramped, lt(n(-kernel.offset), slot)));

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

  /// The stencil of the scheme of axis `d`, an axis of knots, at the coordinate `x`. The knots
  /// below the coordinate (or at it) are counted by a binary search of a fixed number of steps, as
  /// many as a grid of the most knots the axis can have needs, which a NaN leaves at none. Linear
  /// and nearest interpolate in the cell from the last knot at or below the coordinate (a
  /// coordinate on a knot has s = 0 there, and the value of the knot), previous and next take the
  /// value of the last knot below it, or of the first at or above it (a coordinate on a knot
  /// belongs to the cell below it, as in easi's layered model). Beyond the knots the value at the
  /// nearest end.
  Axis locateKnots(GridId grid,
                   std::int32_t d,
                   NodeId x,
                   const df::GridDesc& desc,
                   df::Interpolation interpolation,
                   bool fails) {
    Arena& a = arena();
    if (desc.timeAxis.has_value() && *desc.timeAxis == d) {
      fail("the time axis of the grid " + desc.label() + " has knots; a time axis is uniform");
    }
    if (interpolation == df::Interpolation::Cubic) {
      fail("an axis of knots interpolates nearest, linear, previous or next, not cubic");
    }
    Axis axis;
    const NodeId num = a.gridParam(grid, GridField::Num, d);
    axis.top = isub(num, n(1));
    axis.stride = a.gridParam(grid, GridField::Stride, d);
    const auto knot = [&](NodeId at) { return a.gridLoad(grid, knotsOf(d), at); };
    // the steps of the search: down from the largest power of two not above the most knots
    const auto axisIndex = static_cast<std::size_t>(d);
    const std::size_t most = desc.data != nullptr && axisIndex < desc.data->knots.size()
                                 ? std::max<std::size_t>(1, desc.data->knots[axisIndex].size())
                                 : df::MaxKnotsOfFile;
    std::int64_t highest = 1;
    while (static_cast<std::size_t>(highest) * 2 <= most) {
      highest *= 2;
    }
    const auto counted = [&](bool orAt) {
      NodeId found = n(0);
      for (std::int64_t step = highest; step >= 1; step /= 2) {
        const NodeId candidate = iadd(found, n(step));
        const NodeId probe = knot(isub(imin(candidate, num), n(1)));
        const NodeId taken = orAt ? le(probe, x) : lt(probe, x);
        found = choose(both(le(candidate, num), taken), candidate, found, I64);
      }
      return found;
    };
    if (fails) {
      const NodeId undefined = eq(eq(x, x), truth(false));
      axis.clamped = either(either(lt(x, knot(n(0))), undefined), lt(knot(axis.top), x));
    }

    if (interpolation == df::Interpolation::Previous || interpolation == df::Interpolation::Next) {
      const NodeId below = counted(false);
      axis.start = interpolation == df::Interpolation::Previous ? imax(isub(below, n(1)), n(0))
                                                                : imin(below, axis.top);
      axis.weights = {k(1)};
      return axis;
    }

    // the cell from the last knot at or below the coordinate, the first cell below the first knot
    // and the last one beyond the last knot; one knot is a cell of its own
    const NodeId atOrBelow = counted(true);
    const NodeId index = imax(imin(isub(atOrBelow, n(1)), isub(num, n(2))), n(0));
    const NodeId raw = smul(ssub(x, knot(index)), a.gridLoad(grid, inverseSpacingsOf(d), index));
    // in [0, 1], a NaN at 0, and 1 from the last knot on, where the product need not be one
    const NodeId within = choose(le(k(0), raw), choose(lt(k(1), raw), k(1), raw, F64), k(0), F64);
    const NodeId s = choose(eq(atOrBelow, num), k(1), within, F64);
    if (interpolation == df::Interpolation::Nearest) {
      // half up, as the nearest node of a uniform axis
      axis.start = imin(choose(le(k(0.5), s), iadd(index, n(1)), index, I64), axis.top);
      axis.weights = {k(1)};
      return axis;
    }
    axis.start = index;
    axis.weights = {ssub(k(1), s), s};
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
        fail("the grid " + program.grids()[static_cast<std::size_t>(node.grid)].label() +
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
  case CheckKind::User:
    logError() << "expr: a point does not meet a check of the program on the grid" << path;
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
