// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/Ir.h"
#include "Expr/Program.h"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace seissol::expr {

namespace {

[[noreturn]] void fail(const std::string& what) { throw std::invalid_argument("expr: " + what); }

void checkNames(const Program& program) {
  std::unordered_set<std::string> inputs;
  for (const VarSpec& v : program.inputs()) {
    if (!inputs.insert(v.name).second) {
      fail("duplicate input '" + v.name + "'");
    }
  }
  std::unordered_set<std::string> outputs;
  for (const VarSpec& v : program.outputs()) {
    if (!outputs.insert(v.name).second) {
      fail("duplicate output '" + v.name + "'");
    }
  }
  for (const StateSpec& s : program.state()) {
    // A name that is both an input and a state has two readers with two
    // different storage classes and no rule to pick between them.
    if (inputs.count(s.name) != 0) {
      fail("state '" + s.name + "' collides with an input of the same name");
    }
    if (!inputs.insert(s.name).second) {
      fail("duplicate state '" + s.name + "'");
    }
  }
  if (program.roots().size() != program.outputs().size()) {
    fail("root count (" + std::to_string(program.roots().size()) +
         ") does not match output count (" + std::to_string(program.outputs().size()) + ")");
  }
}

// Every channel a Field may legitimately name, as interned channel ids.
std::unordered_set<int> boundChannels(const Program& program) {
  const Arena& arena = program.arena();
  std::unordered_set<int> bound;
  for (const VarSpec& v : program.inputs()) {
    const int ch = arena.findChannel(v.name);
    if (ch >= 0) {
      bound.insert(ch);
    }
  }
  for (const StateSpec& s : program.state()) {
    const int ch = arena.findChannel(s.name);
    if (ch >= 0) {
      bound.insert(ch);
    }
  }
  return bound;
}

void checkNode(const Program& program, NodeId id, const std::unordered_set<int>& bound) {
  const Arena& arena = program.arena();
  const Node& n = arena[id];

  if (requiresImplicitState(n.kind)) {
    fail(std::string("node kind '") + name(n.kind) +
         "' needs per-point state and a timestep the scripting ABI does not carry; declare a "
         "Program state slot instead");
  }
  if (requiresElementContext(n.kind)) {
    fail(std::string("node kind '") + name(n.kind) +
         "' only means something inside an element / output-cadence context, which a pointwise "
         "program does not have");
  }

  switch (n.kind) {
  case Kind::Const:
    break;
  case Kind::Field:
    if (bound.count(n.ch) == 0) {
      fail("channel '" + arena.channelName(n.ch) +
           "' is read but declared as neither an input "
           "nor a state");
    }
    break;
  case Kind::PW: {
    const int a = arity(n.fn);
    if (n.a == NoNode || (a >= 2 && n.b == NoNode) || (a >= 3 && n.c == NoNode)) {
      fail(std::string("op '") + name(n.fn) + "' has arity " + std::to_string(a) +
           " but not that many operands are set");
    }
    if ((a < 2 && n.b != NoNode) || (a < 3 && n.c != NoNode)) {
      fail(std::string("op '") + name(n.fn) + "' has arity " + std::to_string(a) +
           " but carries extra operands");
    }
    break;
  }
  case Kind::Lookup:
    if (n.grid < 0 || static_cast<std::size_t>(n.grid) >= program.grids().size()) {
      fail("lookup references grid id " + std::to_string(n.grid) + ", but the program declares " +
           std::to_string(program.grids().size()));
    }
    if (n.argCount < 1 || n.argCount > MaxLookupDimension) {
      fail("lookup has " + std::to_string(n.argCount) + " coordinates, expected 1 .. " +
           std::to_string(MaxLookupDimension));
    }
    break;
  case Kind::Dx:
  case Kind::Cumint:
  case Kind::Fold:
  case Kind::Sample:
    break; // already rejected above
  }
}

void walk(const Program& program,
          NodeId root,
          std::vector<bool>& seen,
          const std::unordered_set<int>& bound,
          const char* what,
          std::size_t which) {
  const Arena& arena = program.arena();
  if (root == NoNode || root < 0 || static_cast<std::size_t>(root) >= arena.size()) {
    fail(std::string(what) + " " + std::to_string(which) + " has an out-of-range root id " +
         std::to_string(root));
  }
  // A cycle is unreachable by construction: a node can only reference ids that
  // already existed when it was interned, so the arena is topologically ordered
  // by id. The check that used to be needed on the sderiv side (check_acyclic)
  // belongs to its *pre*-interning surface IR, not here.
  std::vector<NodeId> stack{root};
  std::vector<NodeId> kids;
  while (!stack.empty()) {
    const NodeId id = stack.back();
    stack.pop_back();
    if (seen[id]) {
      continue;
    }
    seen[id] = true;
    checkNode(program, id, bound);
    arena.children(id, kids);
    for (const NodeId k : kids) {
      if (k == NoNode || k < 0 || static_cast<std::size_t>(k) >= arena.size()) {
        fail("node " + std::to_string(id) + " has an out-of-range operand id " + std::to_string(k));
      }
      if (!seen[k]) {
        stack.push_back(k);
      }
    }
  }
}

} // namespace

void validate(const Program& program) {
  checkNames(program);
  const std::unordered_set<int> bound = boundChannels(program);
  std::vector<bool> seen(program.arena().size(), false);

  for (std::size_t i = 0; i < program.roots().size(); ++i) {
    walk(program, program.roots()[i], seen, bound, "output", i);
  }
  for (std::size_t i = 0; i < program.state().size(); ++i) {
    walk(program, program.state()[i].root, seen, bound, "state", i);
  }

  // Grids that nothing reads are not an error — a frontend may intern a grid and
  // then constant-fold the only reference away — but they are worth a note,
  // because the usual cause is a Lookup that silently lost its consumer.
  if (program.roots().empty() && program.state().empty()) {
    fail("program has neither an output nor a state slot");
  }
}

} // namespace seissol::expr
