// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/Rewrite.h"

#include "Expr/Ir.h"
#include "Expr/Program.h"

#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace seissol::expr {

namespace {

/// The nodes the roots and state roots reach. The arena is topologically ordered by id, so a
/// rebuild in ascending id order meets every operand before its users.
std::vector<char> reachable(const Program& program) {
  const Arena& arena = program.arena();
  std::vector<char> seen(arena.size(), 0);
  std::vector<NodeId> stack;
  std::vector<NodeId> kids;
  for (const NodeId root : program.roots()) {
    stack.push_back(root);
  }
  for (const auto& state : program.state()) {
    stack.push_back(state.root);
  }
  while (!stack.empty()) {
    const NodeId id = stack.back();
    stack.pop_back();
    if (id == NoNode || seen[id] != 0) {
      continue;
    }
    seen[id] = 1;
    arena.children(id, kids);
    for (const NodeId kid : kids) {
      stack.push_back(kid);
    }
  }
  return seen;
}

} // namespace

void substituteChannels(Program& program, const std::map<std::string, ChannelBuilder>& builders) {
  const Arena& old = program.arena();
  std::unordered_set<std::string> inputs;
  for (const auto& input : program.inputs()) {
    inputs.insert(input.name);
  }

  Program rebuilt;
  rebuilt.setComputeType(program.computeType());
  for (const auto& grid : program.grids()) {
    rebuilt.internGrid(grid);
  }
  for (const auto& spec : program.matrices()) {
    rebuilt.internMatrix(spec.name, spec.shape);
  }
  for (const auto& spec : program.blocks()) {
    rebuilt.internBlock(spec.name, spec.length);
  }

  // Bottom-up over the reachable nodes, with interning on reinsertion: two reads of the same
  // channel become one contraction node, and so does a contraction an earlier rewrite produced.
  Arena& arena = rebuilt.arena();
  const auto live = reachable(program);
  std::vector<NodeId> map(old.size(), NoNode);
  std::vector<NodeId> coordinates;
  for (NodeId id = 0; id < static_cast<NodeId>(old.size()); ++id) {
    if (live[id] == 0) {
      continue;
    }
    const Node& node = old[id];
    switch (node.kind) {
    case Kind::Const:
      map[id] = arena.konst(node.value);
      break;
    case Kind::Field: {
      const std::string& name = old.channelName(node.ch);
      const auto builder = builders.find(name);
      map[id] = (builder != builders.end() && inputs.count(name) != 0) ? builder->second(arena)
                                                                       : arena.field(name);
      break;
    }
    case Kind::PW:
      switch (arity(node.fn)) {
      case 1:
        map[id] = arena.pw(node.fn, map[node.a]);
        break;
      case 2:
        map[id] = arena.pw(node.fn, map[node.a], map[node.b]);
        break;
      default:
        map[id] = arena.pw(node.fn, map[node.a], map[node.b], map[node.c]);
        break;
      }
      break;
    case Kind::Lookup: {
      coordinates.clear();
      const NodeId* args = old.args(node);
      for (std::int32_t i = 0; i < node.argCount; ++i) {
        coordinates.push_back(map[args[i]]);
      }
      map[id] = arena.lookup(node.grid, node.comp, coordinates);
      break;
    }
    case Kind::Contract:
      map[id] = arena.contract(node.matrix, node.block);
      break;
    case Kind::Dx:
      map[id] = arena.dx(node.axis, map[node.a]);
      break;
    case Kind::Cumint:
      map[id] = arena.cumint(map[node.a]);
      break;
    case Kind::Fold:
      map[id] = arena.fold(node.red, map[node.a]);
      break;
    case Kind::Sample:
      map[id] = arena.sample(map[node.a]);
      break;
    }
  }

  // The signature in its old order, minus the inputs nothing reads any more, plus the channels
  // the builders introduced, in the order they were created.
  std::unordered_set<std::string> states;
  for (const auto& state : program.state()) {
    states.insert(state.name);
  }
  for (const auto& input : program.inputs()) {
    if (arena.findChannel(input.name) >= 0) {
      rebuilt.addInput(input.name, input.type);
    }
  }
  for (std::size_t ch = 0; ch < arena.channelCount(); ++ch) {
    const std::string& name = arena.channelName(static_cast<int>(ch));
    if (inputs.count(name) == 0 && states.count(name) == 0) {
      rebuilt.addInput(name, reader::scripting::DataType::F64);
    }
  }
  for (const auto& state : program.state()) {
    rebuilt.addState(state.name, state.initial, map[state.root]);
  }
  for (std::size_t i = 0; i < program.outputs().size(); ++i) {
    rebuilt.addOutput(
        program.outputs()[i].name, program.outputs()[i].type, map[program.roots()[i]]);
  }

  validate(rebuilt);
  program = std::move(rebuilt);
}

void substituteByContraction(Program& program,
                             MatrixId matrix,
                             const std::map<std::string, BlockId>& blocks) {
  if (matrix < 0 || static_cast<std::size_t>(matrix) >= program.matrices().size()) {
    throw std::invalid_argument("expr: contraction against the undeclared matrix id " +
                                std::to_string(matrix));
  }
  std::map<std::string, ChannelBuilder> builders;
  for (const auto& [name, block] : blocks) {
    builders.emplace(
        name, [matrix, block = block](Arena& arena) { return arena.contract(matrix, block); });
  }
  substituteChannels(program, builders);
}

} // namespace seissol::expr
