// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/Types.h"

#include "Expr/Ir.h"
#include "Expr/Program.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace seissol::expr {

namespace {

[[noreturn]] void fail(const std::string& what) { throw std::invalid_argument("expr: " + what); }

/// The type of a value, or none for a value of literals alone, which takes the type it is used in.
using MaybeType = std::optional<ValueType>;

std::string operationName(const Node& node) {
  if (node.kind == Kind::PW) {
    return std::string("`") + name(node.fn) + "`";
  }
  if (node.kind == Kind::Lookup) {
    return "the lookup of a grid";
  }
  return std::string("`") + name(node.kind) + "`";
}

class Typer {
  public:
  explicit Typer(const Program& program)
      : in_(program), arena_(program.arena()), types_(program.arena().size()),
        reachable_(program.arena().size(), 0) {
    for (const StateSpec& state : program.state()) {
      const int channel = arena_.findChannel(state.name);
      if (channel >= 0) {
        stateTypes_.emplace(channel, state.type);
      }
    }
    markReachable();
    infer();
  }

  /// Whether a value of the program, or a state, is in binary64.
  [[nodiscard]] bool hasF64() const {
    for (std::size_t id = 0; id < types_.size(); ++id) {
      if (reachable_[id] != 0 && types_[id] == ValueType::F64) {
        return true;
      }
    }
    for (const StateSpec& state : in_.state()) {
      if (state.type == ValueType::F64) {
        return true;
      }
    }
    return false;
  }

  [[nodiscard]] Program typedProgram() {
    out_.setComputeType(in_.computeType());
    for (const auto& grid : in_.grids()) {
      out_.internGrid(grid);
    }
    for (const auto& matrix : in_.matrices()) {
      out_.internMatrix(matrix.name, matrix.shape);
    }
    for (const auto& block : in_.blocks()) {
      out_.internBlock(block.name, block.length);
    }
    for (const auto& input : in_.inputs()) {
      out_.addInput(input.name, input.type);
    }
    for (const StateSpec& state : in_.state()) {
      const MaybeType& next = types_[state.root];
      if (next.has_value() && *next != state.type) {
        fail(std::string("the state `") + state.name + "` is " + name(state.type) +
             ", its next value " + name(*next) + "; convert it with " + name(state.type) + "()");
      }
    }
    buildAll();
    for (std::size_t i = 0; i < in_.outputs().size(); ++i) {
      const auto& output = in_.outputs()[i];
      const NodeId root = in_.roots()[i];
      out_.addOutput(
          output.name, output.type, built_.at(keyOf(root, typeIn(root, ValueType::Real))));
    }
    for (const StateSpec& state : in_.state()) {
      out_.addState(state.name,
                    state.initial,
                    built_.at(keyOf(state.root, typeIn(state.root, state.type))),
                    state.type);
    }
    return std::move(out_);
  }

  private:
  void markReachable() {
    std::vector<NodeId> stack;
    std::vector<NodeId> kids;
    for (const NodeId root : in_.roots()) {
      stack.push_back(root);
    }
    for (const StateSpec& state : in_.state()) {
      stack.push_back(state.root);
    }
    while (!stack.empty()) {
      const NodeId id = stack.back();
      stack.pop_back();
      if (id == NoNode || reachable_[id] != 0) {
        continue;
      }
      reachable_[id] = 1;
      arena_.children(id, kids);
      stack.insert(stack.end(), kids.begin(), kids.end());
    }
  }

  static MaybeType unify(const MaybeType& a, const MaybeType& b, const Node& node) {
    if (!a.has_value()) {
      return b;
    }
    if (!b.has_value() || *a == *b) {
      return a;
    }
    fail("the operands of " + operationName(node) + " are " + name(*a) + " and " + name(*b) +
         "; convert one of them with " + name(*a) + "() or " + name(*b) + "()");
  }

  /// The type of every reachable value, from its operands: the arena is in topological order.
  void infer() {
    std::vector<NodeId> kids;
    for (std::size_t index = 0; index < types_.size(); ++index) {
      if (reachable_[index] == 0) {
        continue;
      }
      const auto id = static_cast<NodeId>(index);
      const Node& node = arena_[id];
      switch (node.kind) {
      case Kind::Const:
        types_[index] = std::nullopt;
        break;
      case Kind::Field: {
        // an input in the compute type, unless typed() made it a load in another
        const auto state = stateTypes_.find(node.ch);
        types_[index] = state != stateTypes_.end() ? state->second : node.type;
        break;
      }
      case Kind::Cast:
        types_[index] = node.type;
        break;
      case Kind::Contract:
        types_[index] = ValueType::Real;
        break;
      case Kind::GridParam:
      case Kind::GridLoad:
        // the parts of a written-out lookup have the types they were made with
        types_[index] = node.type;
        break;
      default: {
        // one type for all operands
        arena_.children(id, kids);
        MaybeType type;
        for (const NodeId kid : kids) {
          type = unify(type, types_[kid], node);
        }
        types_[index] = type;
        break;
      }
      }
    }
  }

  [[nodiscard]] static std::uint64_t keyOf(NodeId id, ValueType type) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(id)) << 3U) |
           static_cast<std::uint64_t>(type);
  }

  [[nodiscard]] ValueType typeIn(NodeId id, ValueType context) const {
    return types_[id].value_or(context);
  }

  [[nodiscard]] bool isInput(NodeId id) const {
    const Node& node = arena_[id];
    return node.kind == Kind::Field && stateTypes_.count(node.ch) == 0;
  }

  /// The operands of `id` in type `type`, with the types they are built in.
  void operandsOf(NodeId id, ValueType type, std::vector<std::pair<NodeId, ValueType>>& out) const {
    out.clear();
    const Node& node = arena_[id];
    switch (node.kind) {
    case Kind::Const:
    case Kind::Field:
    case Kind::Contract:
      return;
    case Kind::Cast: {
      const MaybeType& source = types_[node.a];
      if (!source.has_value()) {
        // a literal, in the type converted to
        out.emplace_back(node.a, node.type);
      } else if (*source == node.type || !isInput(node.a)) {
        out.emplace_back(node.a, *source);
      }
      // an input converted to another type is a load in that type, without operands
      return;
    }
    default: {
      std::vector<NodeId> kids;
      arena_.children(id, kids);
      for (const NodeId kid : kids) {
        out.emplace_back(kid, type);
      }
      return;
    }
    }
  }

  /// The node of `id` in type `type`, of the operands built before.
  [[nodiscard]] NodeId make(NodeId id, ValueType type) {
    const Node& node = arena_[id];
    Arena& arena = out_.arena();
    const auto operand = [&](NodeId kid, ValueType kidType) {
      return built_.at(keyOf(kid, kidType));
    };
    switch (node.kind) {
    case Kind::Const:
      return arena.konst(node.value, type);
    case Kind::Field:
      return arena.field(arena.channel(arena_.channelName(node.ch)), type);
    case Kind::Cast: {
      const MaybeType& source = types_[node.a];
      if (!source.has_value()) {
        return operand(node.a, node.type);
      }
      if (*source == node.type) {
        return operand(node.a, *source);
      }
      if (isInput(node.a)) {
        return arena.field(arena.channel(arena_.channelName(arena_[node.a].ch)), node.type);
      }
      return arena.cast(node.type, operand(node.a, *source));
    }
    case Kind::PW: {
      const int count = arity(node.fn);
      if (node.strict) {
        return arena.strictPw(node.fn, operand(node.a, type), operand(node.b, type), type);
      }
      if (count == 1) {
        return arena.pw(node.fn, operand(node.a, type), type);
      }
      if (count == 2) {
        return arena.pw(node.fn, operand(node.a, type), operand(node.b, type), type);
      }
      return arena.pw(
          node.fn, operand(node.a, type), operand(node.b, type), operand(node.c, type), type);
    }
    case Kind::Lookup: {
      std::vector<NodeId> coordinates;
      const NodeId* args = arena_.args(node);
      for (std::int32_t i = 0; i < node.argCount; ++i) {
        coordinates.push_back(operand(args[i], type));
      }
      return arena.lookup(node.grid, node.comp, coordinates, type);
    }
    case Kind::Contract:
      return arena.contract(node.matrix, node.block, ValueType::Real);
    case Kind::GridParam:
      return arena.gridParam(node.grid, static_cast<GridField>(node.comp), node.axis);
    case Kind::GridLoad:
      return arena.gridLoad(node.grid, node.comp, operand(node.a, ValueType::I64));
    }
    fail(std::string("node kind '") + name(node.kind) + "' has no type");
  }

  /// The nodes of the values the outputs and the states need, each in the types it is needed in
  /// (a literal may be needed in several): found from the roots down, then made in the order of the
  /// arena of the program, so that a program without conversions keeps the order of its values,
  /// and the lowering emits its instructions as before (a device kernel that contracts into fused
  /// multiply-adds then keeps its bits).
  void buildAll() {
    const std::size_t count = types_.size();
    // per node, a bit per type it is needed in
    std::vector<std::uint8_t> needed(count, 0);
    const auto bit = [](ValueType type) {
      return static_cast<std::uint8_t>(1U << static_cast<unsigned>(type));
    };
    std::vector<std::pair<NodeId, ValueType>> stack;
    for (const NodeId root : in_.roots()) {
      stack.emplace_back(root, typeIn(root, ValueType::Real));
    }
    for (const StateSpec& state : in_.state()) {
      stack.emplace_back(state.root, typeIn(state.root, state.type));
    }
    std::vector<std::pair<NodeId, ValueType>> operands;
    while (!stack.empty()) {
      const auto [id, type] = stack.back();
      stack.pop_back();
      auto& mark = needed[static_cast<std::size_t>(id)];
      if ((mark & bit(type)) != 0) {
        continue;
      }
      mark = static_cast<std::uint8_t>(mark | bit(type));
      operandsOf(id, type, operands);
      stack.insert(stack.end(), operands.begin(), operands.end());
    }
    for (std::size_t index = 0; index < count; ++index) {
      for (const auto type : {ValueType::Real, ValueType::F32, ValueType::F64}) {
        if ((needed[index] & bit(type)) != 0) {
          const auto id = static_cast<NodeId>(index);
          built_.emplace(keyOf(id, type), make(id, type));
        }
      }
    }
  }

  const Program& in_;
  const Arena& arena_;
  std::vector<MaybeType> types_;
  std::vector<char> reachable_;
  std::unordered_map<int, ValueType> stateTypes_;
  Program out_;
  std::unordered_map<std::uint64_t, NodeId> built_;
};

} // namespace

ComputeType resolve(ValueType type, ComputeType computeType) {
  switch (type) {
  case ValueType::Real:
    return computeType;
  case ValueType::F32:
    return ComputeType::F32;
  case ValueType::F64:
    return ComputeType::F64;
  case ValueType::I64:
    return ComputeType::I64;
  case ValueType::Bool:
    return ComputeType::Bool;
  }
  return computeType;
}

Program typed(const Program& program) { return Typer(program).typedProgram(); }

ComputeType dataType(const Program& program) {
  if (program.computeType() == ComputeType::F64) {
    return ComputeType::F64;
  }
  return Typer(program).hasF64() ? ComputeType::F64 : ComputeType::F32;
}

} // namespace seissol::expr
