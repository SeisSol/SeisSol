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

[[noreturn]] void fail(const std::string& what) { throw std::invalid_argument("expr: " + what); }

/// One axis of a lookup, located: the first node of the stencil, the weights of its nodes, the last
/// node of the axis and the elements between two nodes; and whether the coordinate was clamped.
struct Axis {
  NodeId start{NoNode};
  std::vector<NodeId> weights;
  NodeId top{NoNode};
  NodeId stride{NoNode};
  NodeId clamped{NoNode};
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

  // the whole numbers, which are exact in binary64: their order and contraction do not matter
  NodeId k(double value) { return arena().konst(value, F64); }
  NodeId add(NodeId x, NodeId y) { return arena().pw(Fn::Add, x, y, F64); }
  NodeId sub(NodeId x, NodeId y) { return arena().pw(Fn::Sub, x, y, F64); }
  NodeId mul(NodeId x, NodeId y) { return arena().pw(Fn::Mul, x, y, F64); }
  // the arithmetic of the samples, each operation rounded on its own
  NodeId sadd(NodeId x, NodeId y) { return arena().strictPw(Fn::Add, x, y, F64); }
  NodeId ssub(NodeId x, NodeId y) { return arena().strictPw(Fn::Sub, x, y, F64); }
  NodeId smul(NodeId x, NodeId y) { return arena().strictPw(Fn::Mul, x, y, F64); }
  // comparisons, which give 1 or 0, and the choices between two values
  NodeId lt(NodeId x, NodeId y) { return arena().pw(Fn::Lt, x, y, F64); }
  NodeId le(NodeId x, NodeId y) { return arena().pw(Fn::Le, x, y, F64); }
  NodeId eq(NodeId x, NodeId y) { return arena().pw(Fn::Eq, x, y, F64); }
  NodeId both(NodeId x, NodeId y) { return arena().pw(Fn::And, x, y, F64); }
  NodeId either(NodeId x, NodeId y) { return arena().pw(Fn::Or, x, y, F64); }
  NodeId choose(NodeId condition, NodeId x, NodeId y) {
    return arena().pw(Fn::Select, condition, x, y, F64);
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
      return a.gridLoad(node.grid, mapped(node.a));
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
      axes.push_back(locate(grid, d, wide ? coordinate : a.cast(F64, coordinate), desc, fails));
    }

    const NodeId base = mul(k(node.comp), a.gridParam(grid, GridField::ComponentStride, 0));
    const NodeId sample = reduce(grid, axes, 0, base);

    if (fails) {
      NodeId outside = axes.front().clamped;
      for (std::size_t d = 1; d < axes.size(); ++d) {
        outside = either(outside, axes[d].clamped);
      }
      result_.checks.push_back(GridCheck{grid, outside});
    }
    return wide ? sample : a.cast(node.type, sample);
  }

  /// The stencil of the scheme on axis `d` at the coordinate `x`: locateAxis in Interpolation.cpp.
  Axis locate(GridId grid, std::int32_t d, NodeId x, const df::GridDesc& desc, bool fails) {
    Arena& a = arena();
    Axis axis;
    const NodeId num = a.gridParam(grid, GridField::Num, d);
    axis.top = sub(num, k(1));
    axis.stride = a.gridParam(grid, GridField::Stride, d);
    const NodeId raw = smul(ssub(x, a.gridParam(grid, GridField::Min, d)),
                            a.gridParam(grid, GridField::DeltaInv, d));
    // `raw > 0`, not `raw >= 0`: NaN clamps onto the lower end
    const NodeId above = lt(axis.top, raw);
    const NodeId clampedRaw = choose(lt(k(0), raw), choose(above, axis.top, raw), k(0));
    if (fails) {
      axis.clamped = either(either(lt(raw, k(0)), eq(eq(raw, raw), k(0))), above);
    }

    if (desc.interpolation == df::Interpolation::Nearest) {
      // half away from zero, as std::lround
      axis.start = a.pw(Fn::Round, clampedRaw, F64);
      axis.weights = {k(1)};
      return axis;
    }

    // the base cell is the last one with a right neighbour: a coordinate on the last node is s = 1
    // in the cell before
    const NodeId twoNodes = le(k(2), num);
    const NodeId maxBase = choose(twoNodes, sub(num, k(2)), k(0));
    const NodeId lower = a.pw(Fn::Floor, clampedRaw, F64);
    const NodeId past = lt(maxBase, lower);
    const NodeId index = choose(past, maxBase, lower);
    const NodeId s = choose(past, choose(twoNodes, k(1), k(0)), ssub(clampedRaw, lower));
    const NodeId oneMinusS = ssub(k(1), s);

    if (desc.interpolation == df::Interpolation::Linear) {
      axis.start = index;
      axis.weights = {oneMinusS, s};
      return axis;
    }

    // cubic: four nodes from the one before the base cell, moved inside the axis at its ends,
    // where the weights fall back to linear ones between the base and its neighbour
    const df::StencilKernel kernel = df::kernelOf(desc.interpolation);
    const auto width = static_cast<double>(kernel.width);
    const NodeId fourNodes = le(k(width), num);
    const NodeId maxStart = sub(num, k(width));
    const NodeId first = choose(fourNodes, add(index, k(static_cast<double>(kernel.offset))), k(0));
    const NodeId start =
        choose(lt(first, k(0)),
               k(0),
               choose(lt(maxStart, first), choose(lt(k(0), maxStart), maxStart, k(0)), first));
    const NodeId slot = sub(index, start);
    const NodeId full = both(fourNodes, eq(slot, k(static_cast<double>(-kernel.offset))));

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
      const NodeId fallback = choose(eq(slot, k(static_cast<double>(j))),
                                     oneMinusS,
                                     choose(eq(slot, k(static_cast<double>(j) - 1)), s, k(0)));
      axis.weights.push_back(choose(full, cubic[j], fallback));
    }
    return axis;
  }

  /// The tensor product over the axes from `d` on, of the elements from `offset` on: the sum over
  /// the last axis innermost, each sum from zero in the order of its nodes, as the reduction of
  /// sampleBatch() forms them.
  NodeId reduce(GridId grid, const std::vector<Axis>& axes, std::size_t d, NodeId offset) {
    Arena& a = arena();
    const Axis& axis = axes[d];
    NodeId sum = k(0);
    for (std::size_t i = 0; i < axis.weights.size(); ++i) {
      // only the upper end needs a clamp: the start is not negative
      const NodeId shifted = i == 0 ? axis.start : add(axis.start, k(static_cast<double>(i)));
      const NodeId node = a.pw(Fn::Min, shifted, axis.top, F64);
      const NodeId at = add(offset, mul(node, axis.stride));
      const NodeId value =
          d + 1 == axes.size() ? a.gridLoad(grid, at) : reduce(grid, axes, d + 1, at);
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
  }
  return 0.0;
}

} // namespace seissol::expr
