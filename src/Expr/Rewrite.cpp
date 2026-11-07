// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/Rewrite.h"

#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Reader/Scripting/DataTable.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
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

NodeId replaceChannels(Arena& arena,
                       NodeId root,
                       const std::unordered_map<int, NodeId>& replacements,
                       std::unordered_map<NodeId, NodeId>& memo) {
  if (root == NoNode) {
    return root;
  }
  const auto found = memo.find(root);
  if (found != memo.end()) {
    return found->second;
  }
  // a copy: interning a node below may move the nodes of the arena
  const Node node = arena[root];
  NodeId result = root;
  switch (node.kind) {
  case Kind::Const:
  case Kind::Contract:
    break;
  case Kind::Field: {
    const auto replacement = replacements.find(node.ch);
    if (replacement != replacements.end()) {
      result = replacement->second;
    }
    break;
  }
  case Kind::PW: {
    const int count = arity(node.fn);
    const NodeId a = replaceChannels(arena, node.a, replacements, memo);
    const NodeId b = count > 1 ? replaceChannels(arena, node.b, replacements, memo) : NoNode;
    const NodeId c = count > 2 ? replaceChannels(arena, node.c, replacements, memo) : NoNode;
    if (a != node.a || b != node.b || c != node.c) {
      if (count == 1) {
        result = arena.pw(node.fn, a);
      } else if (count == 2) {
        result = arena.pw(node.fn, a, b);
      } else {
        result = arena.pw(node.fn, a, b, c);
      }
    }
    break;
  }
  case Kind::Lookup: {
    std::vector<NodeId> coordinates(arena.args(node), arena.args(node) + node.argCount);
    bool changed = false;
    for (auto& coordinate : coordinates) {
      const NodeId replaced = replaceChannels(arena, coordinate, replacements, memo);
      changed |= replaced != coordinate;
      coordinate = replaced;
    }
    if (changed) {
      result = arena.lookup(node.grid, node.comp, coordinates);
    }
    break;
  }
  default:
    throw std::invalid_argument(std::string("cannot replace the channels below a node of kind ") +
                                name(node.kind));
  }
  memo.emplace(root, result);
  return result;
}

namespace {

/// The grids, matrices and blocks of a source program, as they are declared in the program its
/// nodes are added to: declared there when first used, unless declared beforehand.
struct Imports {
  explicit Imports(const Program& source) : source(&source) {}
  const Program* source;
  std::unordered_map<GridId, GridId> grids;
  std::unordered_map<MatrixId, MatrixId> matrices;
  std::unordered_map<BlockId, BlockId> blocks;
  std::unordered_map<NodeId, NodeId> nodes;

  GridId grid(Program& target, GridId id) {
    const auto found = grids.find(id);
    if (found != grids.end()) {
      return found->second;
    }
    const auto imported = target.internGrid(source->grids().at(static_cast<std::size_t>(id)));
    grids.emplace(id, imported);
    return imported;
  }
  MatrixId matrix(Program& target, MatrixId id) {
    const auto found = matrices.find(id);
    if (found != matrices.end()) {
      return found->second;
    }
    const auto& spec = source->matrices().at(static_cast<std::size_t>(id));
    const auto imported = target.internMatrix(spec.name, spec.shape);
    matrices.emplace(id, imported);
    return imported;
  }
  BlockId block(Program& target, BlockId id) {
    const auto found = blocks.find(id);
    if (found != blocks.end()) {
      return found->second;
    }
    const auto& spec = source->blocks().at(static_cast<std::size_t>(id));
    const auto imported = target.internBlock(spec.name, spec.length);
    blocks.emplace(id, imported);
    return imported;
  }
};

/// `root` of the program `imports.source`, built in `target`, whose channels are those of the same
/// names.
NodeId importNode(Program& target, NodeId root, Imports& imports) {
  if (root == NoNode) {
    return root;
  }
  const auto found = imports.nodes.find(root);
  if (found != imports.nodes.end()) {
    return found->second;
  }
  const Arena& source = imports.source->arena();
  const Node node = source[root];
  NodeId result = NoNode;
  switch (node.kind) {
  case Kind::Const:
    result = target.arena().konst(node.value);
    break;
  case Kind::Field:
    result = target.arena().field(source.channelName(node.ch));
    break;
  case Kind::PW: {
    const int count = arity(node.fn);
    const NodeId a = importNode(target, node.a, imports);
    const NodeId b = count > 1 ? importNode(target, node.b, imports) : NoNode;
    const NodeId c = count > 2 ? importNode(target, node.c, imports) : NoNode;
    if (count == 1) {
      result = target.arena().pw(node.fn, a);
    } else if (count == 2) {
      result = target.arena().pw(node.fn, a, b);
    } else {
      result = target.arena().pw(node.fn, a, b, c);
    }
    break;
  }
  case Kind::Lookup: {
    std::vector<NodeId> coordinates(source.args(node), source.args(node) + node.argCount);
    for (auto& coordinate : coordinates) {
      coordinate = importNode(target, coordinate, imports);
    }
    const auto grid = imports.grid(target, node.grid);
    result = target.arena().lookup(grid, node.comp, coordinates);
    break;
  }
  case Kind::Contract: {
    const auto matrix = imports.matrix(target, node.matrix);
    const auto block = imports.block(target, node.block);
    result = target.arena().contract(matrix, block);
    break;
  }
  default:
    throw std::invalid_argument(std::string("cannot rebuild a program with a node of kind ") +
                                name(node.kind));
  }
  imports.nodes.emplace(root, result);
  return result;
}

template <typename Spec>
const Spec* findByName(const std::vector<Spec>& specs, const std::string& name) {
  for (const auto& spec : specs) {
    if (spec.name == name) {
      return &spec;
    }
  }
  return nullptr;
}

} // namespace

Program mergePrograms(const Program& first, const Program& second) {
  for (const auto& output : second.outputs()) {
    if (findByName(first.outputs(), output.name) != nullptr) {
      throw std::invalid_argument("both give the output `" + output.name + "`");
    }
  }
  for (const auto& state : second.state()) {
    if (findByName(first.state(), state.name) != nullptr) {
      throw std::invalid_argument("both declare the state `" + state.name + "`");
    }
    if (findByName(first.inputs(), state.name) != nullptr) {
      throw std::invalid_argument("`" + state.name +
                                  "` is a state of one and an input of the other");
    }
  }
  for (const auto& state : first.state()) {
    if (findByName(second.inputs(), state.name) != nullptr) {
      throw std::invalid_argument("`" + state.name +
                                  "` is a state of one and an input of the other");
    }
  }
  for (const auto& input : second.inputs()) {
    const auto* other = findByName(first.inputs(), input.name);
    if (other != nullptr && other->type != input.type) {
      throw std::invalid_argument("both read `" + input.name + "`, with different types");
    }
  }

  Program merged = first;
  if (second.computeType() == ComputeType::F64) {
    merged.setComputeType(ComputeType::F64);
  }
  Imports imports(second);
  for (std::size_t i = 0; i < second.grids().size(); ++i) {
    imports.grid(merged, static_cast<GridId>(i));
  }
  for (std::size_t i = 0; i < second.matrices().size(); ++i) {
    imports.matrix(merged, static_cast<MatrixId>(i));
  }
  for (std::size_t i = 0; i < second.blocks().size(); ++i) {
    imports.block(merged, static_cast<BlockId>(i));
  }
  for (const auto& input : second.inputs()) {
    if (findByName(merged.inputs(), input.name) == nullptr) {
      merged.addInput(input.name, input.type);
    }
  }
  for (std::size_t i = 0; i < second.outputs().size(); ++i) {
    const auto& output = second.outputs()[i];
    merged.addOutput(output.name, output.type, importNode(merged, second.roots()[i], imports));
  }
  for (const auto& state : second.state()) {
    merged.addState(state.name, state.initial, importNode(merged, state.root, imports));
  }
  validate(merged);
  return merged;
}

Program restrictProgram(const Program& program,
                        const std::function<bool(const std::string&)>& available,
                        const std::function<bool(const std::string&)>& wanted,
                        Restriction* restriction) {
  const auto& states = program.state();
  const auto isState = [&](const std::string& name) { return findByName(states, name) != nullptr; };

  // what the update of each state and each output reads
  std::vector<std::set<std::string>> stateReads;
  stateReads.reserve(states.size());
  for (const auto& state : states) {
    stateReads.push_back(channelsRead(program, {state.root}));
  }
  std::vector<std::set<std::string>> outputReads;
  outputReads.reserve(program.outputs().size());
  for (const auto root : program.roots()) {
    outputReads.push_back(channelsRead(program, {root}));
  }

  // the states that can be computed: the ones that read only inputs that are there and states
  // that can be computed; and why the others cannot
  std::map<std::string, std::string> missing;
  bool changed = true;
  while (changed) {
    changed = false;
    for (std::size_t s = 0; s < states.size(); ++s) {
      if (missing.count(states[s].name) > 0) {
        continue;
      }
      for (const auto& read : stateReads[s]) {
        const auto lacking = isState(read) ? missing.find(read) : missing.end();
        if (isState(read) ? lacking != missing.end() : !available(read)) {
          missing.emplace(states[s].name, isState(read) ? lacking->second : read);
          changed = true;
          break;
        }
      }
    }
  }
  const auto lackingOf = [&](const std::set<std::string>& reads) -> std::string {
    for (const auto& read : reads) {
      if (isState(read)) {
        const auto lacking = missing.find(read);
        if (lacking != missing.end()) {
          return lacking->second;
        }
      } else if (!available(read)) {
        return read;
      }
    }
    return {};
  };

  // the outputs kept, and the states they need
  std::vector<std::size_t> kept;
  std::set<std::string> needed;
  std::vector<std::string> pending;
  for (std::size_t i = 0; i < program.outputs().size(); ++i) {
    const auto& name = program.outputs()[i].name;
    if (!wanted(name)) {
      if (restriction != nullptr) {
        restriction->unwanted.push_back(name);
      }
      continue;
    }
    const auto lacking = lackingOf(outputReads[i]);
    if (!lacking.empty()) {
      if (restriction != nullptr) {
        restriction->unavailable.emplace_back(name, lacking);
      }
      continue;
    }
    kept.push_back(i);
    for (const auto& read : outputReads[i]) {
      if (isState(read) && needed.insert(read).second) {
        pending.push_back(read);
      }
    }
  }
  while (!pending.empty()) {
    const auto name = pending.back();
    pending.pop_back();
    for (std::size_t s = 0; s < states.size(); ++s) {
      if (states[s].name != name) {
        continue;
      }
      for (const auto& read : stateReads[s]) {
        if (isState(read) && needed.insert(read).second) {
          pending.push_back(read);
        }
      }
    }
  }

  // the inputs the rest reads, in their order
  std::set<std::string> read;
  for (const auto i : kept) {
    read.insert(outputReads[i].begin(), outputReads[i].end());
  }
  for (std::size_t s = 0; s < states.size(); ++s) {
    if (needed.count(states[s].name) > 0) {
      read.insert(stateReads[s].begin(), stateReads[s].end());
    }
  }

  Program restricted;
  restricted.setComputeType(program.computeType());
  for (const auto& input : program.inputs()) {
    if (read.count(input.name) > 0) {
      restricted.addInput(input.name, input.type);
    }
  }
  Imports imports(program);
  for (const auto i : kept) {
    const auto& output = program.outputs()[i];
    restricted.addOutput(
        output.name, output.type, importNode(restricted, program.roots()[i], imports));
  }
  for (const auto& state : states) {
    if (needed.count(state.name) > 0) {
      restricted.addState(state.name, state.initial, importNode(restricted, state.root, imports));
    }
  }
  // nothing left is a program that computes nothing, which its consumer skips
  if (!restricted.outputs().empty() || !restricted.state().empty()) {
    validate(restricted);
  }
  return restricted;
}

} // namespace seissol::expr
