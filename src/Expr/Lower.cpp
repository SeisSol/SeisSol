// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/Lower.h"

#include "Expr/Expand.h"
#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Expr/Types.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace seissol::expr {

namespace {

constexpr std::uint64_t FnvOffset = 0xCBF29CE484222325ULL;
constexpr std::uint64_t FnvPrime = 0x100000001B3ULL;
constexpr std::int64_t CostCap = 1 << 20; // saturate; a diamond DAG is exponential otherwise
constexpr std::int32_t NoSlot = -1;

[[noreturn]] void fail(const std::string& what) { throw std::invalid_argument("expr: " + what); }

// Rough relative cost of recomputing one node. The absolute scale is meaningless
// and only the ratios against LowerOptions::hoistThreshold matter. An element of a
// grid is weighted as a random access, which a written-out lookup has width^d of,
// the most expensive thing this IR can do per point. A contraction costs one
// multiply-add per coefficient.
std::int64_t weight(const Program& program, const Node& n) {
  switch (n.kind) {
  case Kind::Const:
    return 0;
  case Kind::Field:
  case Kind::GridParam:
    return 1;
  case Kind::GridLoad:
    return 4;
  case Kind::Contract:
    return static_cast<std::int64_t>(program.matrices()[n.matrix].shape.cols);
  case Kind::PW:
    switch (n.fn) {
    case Fn::Exp:
    case Fn::Log:
    case Fn::Log2:
    case Fn::Log10:
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
    case Fn::Pow:
      return 8;
    case Fn::Div:
    case Fn::Rcp:
    case Fn::Mod:
    case Fn::Atan2:
    case Fn::Sqrt:
      return 4;
    default:
      return 1;
    }
  default:
    return 1;
  }
}

Space join(Space x, Space y, NodeId where) {
  if (x == Space::Uniform) {
    return y;
  }
  if (y == Space::Uniform || x == y) {
    return x;
  }
  fail("node " + std::to_string(where) + " mixes index spaces (" + name(x) + " and " + name(y) +
       "); an explicit transform is needed between them");
}

// Maps a Field's interned channel id onto either an input slot or a state slot.
struct ChannelMap {
  std::unordered_map<int, std::int32_t> input;
  std::unordered_map<int, std::int32_t> state;

  explicit ChannelMap(const Program& program) {
    const Arena& arena = program.arena();
    for (std::size_t i = 0; i < program.inputs().size(); ++i) {
      const int ch = arena.findChannel(program.inputs()[i].name);
      if (ch >= 0) {
        input.emplace(ch, static_cast<std::int32_t>(i));
      }
    }
    for (std::size_t i = 0; i < program.state().size(); ++i) {
      const int ch = arena.findChannel(program.state()[i].name);
      if (ch >= 0) {
        state.emplace(ch, static_cast<std::int32_t>(i));
      }
    }
  }
};

// Per-node analysis results, indexed by NodeId. Arena ids are already a valid
// topological order — a node can only reference ids that existed when it was
// interned — so every pass here is a single ascending sweep and needs no sort.
struct Analysis {
  std::vector<char> reachable;
  std::vector<char> invariant;
  std::vector<std::int64_t> cost;
  std::vector<Space> space;
};

Analysis analyse(const Program& program,
                 const LowerOptions& options,
                 const ChannelMap& channels,
                 const std::vector<GridCheck>& checks) {
  const Arena& arena = program.arena();
  const std::size_t n = arena.size();
  Analysis an;
  an.reachable.assign(n, 0);
  an.invariant.assign(n, 0);
  an.cost.assign(n, 0);
  an.space.assign(n, Space::Point);

  std::unordered_set<std::string> invariantInputs(options.invariantInputs.begin(),
                                                  options.invariantInputs.end());
  std::unordered_set<GridId> invariantGrids(options.invariantGrids.begin(),
                                            options.invariantGrids.end());

  std::vector<NodeId> stack;
  std::vector<NodeId> kids;
  const auto push = [&](NodeId id) {
    if (id != NoNode && an.reachable[id] == 0) {
      stack.push_back(id);
    }
  };
  for (const NodeId r : program.roots()) {
    push(r);
  }
  for (const StateSpec& s : program.state()) {
    push(s.root);
  }
  for (const GridCheck& check : checks) {
    push(check.root);
  }
  while (!stack.empty()) {
    const NodeId id = stack.back();
    stack.pop_back();
    if (an.reachable[id] != 0) {
      continue;
    }
    an.reachable[id] = 1;
    arena.children(id, kids);
    for (const NodeId k : kids) {
      push(k);
    }
  }

  for (NodeId id = 0; id < static_cast<NodeId>(n); ++id) {
    if (an.reachable[id] == 0) {
      continue;
    }
    const Node& node = arena[id];
    arena.children(id, kids);

    bool inv = true;
    Space sp = Space::Uniform;
    std::int64_t c = weight(program, node);
    for (const NodeId k : kids) {
      inv = inv && an.invariant[k] != 0;
      sp = join(sp, an.space[k], id);
      c = std::min(CostCap, c + an.cost[k]);
    }

    switch (node.kind) {
    case Kind::Const:
      inv = true;
      sp = Space::Uniform;
      break;
    case Kind::Field: {
      sp = Space::Point;
      if (channels.state.count(node.ch) != 0) {
        inv = false; // a state slot changes by definition
      } else {
        const auto it = channels.input.find(node.ch);
        if (it == channels.input.end()) {
          fail("channel '" + arena.channelName(node.ch) +
               "' is bound to neither an input nor a "
               "state; run validate() first");
        }
        inv = invariantInputs.count(program.inputs()[it->second].name) != 0;
      }
      break;
    }
    case Kind::GridParam:
      // the geometry of a grid the store does not update stays
      inv = invariantGrids.count(node.grid) != 0;
      sp = Space::Uniform;
      break;
    case Kind::GridLoad:
      inv = inv && invariantGrids.count(node.grid) != 0;
      break;
    case Kind::Cast:
      break;
    case Kind::Contract:
      // the blocks are the evolving state of the solver, never invariant across calls
      inv = false;
      sp = Space::Point;
      break;
    case Kind::PW:
      break;
    default:
      fail(std::string("node kind '") + name(node.kind) +
           "' cannot be lowered for the scripting path; run validate() first");
    }

    an.invariant[id] = inv ? 1 : 0;
    an.space[id] = sp;
    an.cost[id] = c;
  }
  return an;
}

// Straight-line liveness. Everything is in one basic block per stage, so a
// value's live range is [definition, last use] and slot assignment is a linear
// scan over a free list — no interval tree, no interference graph.
// A slot holds values of one type only; the free slots are kept per type.
class SlotAllocator {
  public:
  explicit SlotAllocator(std::size_t arenaSize) : slotOf_(arenaSize, NoSlot) {}

  [[nodiscard]] std::int32_t slotOf(NodeId id) const { return slotOf_[id]; }

  void release(NodeId id) {
    const std::int32_t s = slotOf_[id];
    if (s != NoSlot) {
      free_.at(index(types_[s])).push_back(s);
      slotOf_[id] = NoSlot;
    }
  }

  std::int32_t allocate(NodeId id, ComputeType type) {
    auto& free = free_.at(index(type));
    std::int32_t s = 0;
    if (free.empty()) {
      s = static_cast<std::int32_t>(types_.size());
      types_.push_back(type);
    } else {
      s = free.back();
      free.pop_back();
    }
    slotOf_[id] = s;
    return s;
  }

  [[nodiscard]] std::int32_t high() const { return static_cast<std::int32_t>(types_.size()); }
  [[nodiscard]] const std::vector<ComputeType>& types() const { return types_; }

  private:
  static std::size_t index(ComputeType type) { return static_cast<std::size_t>(type); }

  std::vector<std::int32_t> slotOf_;
  std::array<std::vector<std::int32_t>, 4> free_;
  std::vector<ComputeType> types_;
};

class Emitter {
  public:
  Emitter(const Program& program,
          const Analysis& analysis,
          const ChannelMap& channels,
          std::vector<std::int32_t>& operands)
      : program_(program), analysis_(analysis), channels_(channels), operands_(operands) {}

  // `values` must be in ascending NodeId order (a valid topological order) and
  // closed under operands, except at the imported-value boundary.
  void emit(const std::vector<NodeId>& values,
            const std::vector<NodeId>& liveOut,
            StageCode& stage) {
    const Arena& arena = program_.arena();
    SlotAllocator alloc(arena.size());

    // Last use of each value inside this stage. Imported values contribute no
    // uses — their operands belong to the stage that produced them. Anything in
    // liveOut survives to the stores at the end and must not be recycled before.
    std::unordered_map<NodeId, std::size_t> lastUse;
    std::vector<NodeId> kids;
    for (std::size_t i = 0; i < values.size(); ++i) {
      const NodeId id = values[i];
      if (isImported(id)) {
        continue;
      }
      arena.children(id, kids);
      for (const NodeId k : kids) {
        lastUse[k] = i;
      }
    }
    for (const NodeId id : liveOut) {
      lastUse[id] = std::numeric_limits<std::size_t>::max();
    }

    for (std::size_t i = 0; i < values.size(); ++i) {
      const NodeId id = values[i];
      const Node& node = arena[id];

      Instruction inst;
      inst.node = id;
      inst.space = analysis_.space[id];
      inst.type = resolve(node.type, program_.computeType());

      const std::int32_t boundary = boundarySlot(id);
      if (boundary >= 0) {
        // Value produced by an earlier stage: bring it in with an explicit load.
        inst.op = Opcode::LoadPersistent;
        inst.imm = boundary;
      } else {
        switch (node.kind) {
        case Kind::Const:
          inst.op = Opcode::Const;
          inst.value = node.value;
          break;
        case Kind::Field: {
          const auto st = channels_.state.find(node.ch);
          if (st != channels_.state.end()) {
            inst.op = Opcode::LoadPersistent;
            inst.imm = st->second; // state slots occupy the first persistent slots
          } else {
            inst.op = Opcode::LoadInput;
            inst.imm = channels_.input.at(node.ch);
          }
          break;
        }
        case Kind::PW:
          inst.op = Opcode::Pw;
          inst.fn = node.fn;
          inst.strict = node.strict;
          break;
        case Kind::GridParam:
          inst.op = Opcode::GridParam;
          inst.grid = node.grid;
          inst.comp = node.comp;
          inst.imm = node.axis;
          break;
        case Kind::GridLoad:
          inst.op = Opcode::GridLoad;
          inst.grid = node.grid;
          inst.comp = node.comp;
          break;
        case Kind::Contract:
          inst.op = Opcode::Contract;
          inst.matrix = node.matrix;
          inst.block = node.block;
          break;
        case Kind::Cast:
          inst.op = Opcode::Convert;
          break;
        default:
          fail(std::string("node kind '") + name(node.kind) + "' reached the emitter");
        }
      }

      if (inst.op == Opcode::Pw || inst.op == Opcode::GridLoad || inst.op == Opcode::Convert) {
        arena.children(id, kids);
        inst.operandBegin = static_cast<std::int32_t>(operands_.size());
        inst.operandCount = static_cast<std::int32_t>(kids.size());
        for (const NodeId k : kids) {
          const std::int32_t s = alloc.slotOf(k);
          if (s == NoSlot) {
            fail("operand " + std::to_string(k) + " of node " + std::to_string(id) +
                 " has no live slot; the stage value set is not closed");
          }
          operands_.push_back(s);
        }
      }

      // Freeing dead operands before allocating the destination lets an op write
      // over an operand it just consumed: every instruction is elementwise -- one
      // lane reads and writes the same index -- now that a lookup is written out.
      const bool inPlaceSafe = true;
      if (inPlaceSafe) {
        releaseDead(alloc, id, lastUse, i);
      }
      inst.dst = alloc.allocate(id, inst.type);
      if (!inPlaceSafe) {
        releaseDead(alloc, id, lastUse, i);
      }

      stage.code.push_back(inst);
    }

    stage.slotCount = alloc.high();
    stage.slotTypes = alloc.types();
    slotOfLiveOut_.clear();
    for (const NodeId id : liveOut) {
      slotOfLiveOut_.emplace(id, alloc.slotOf(id));
    }
  }

  [[nodiscard]] std::int32_t liveOutSlot(NodeId id) const { return slotOfLiveOut_.at(id); }

  // Values that this stage imports rather than computes.
  void setBoundary(const std::unordered_map<NodeId, std::int32_t>& boundary) {
    boundary_ = boundary;
  }

  private:
  [[nodiscard]] std::int32_t boundarySlot(NodeId id) const {
    const auto it = boundary_.find(id);
    return it == boundary_.end() ? -1 : it->second;
  }
  [[nodiscard]] bool isImported(NodeId id) const { return boundary_.count(id) != 0; }

  void releaseDead(SlotAllocator& alloc,
                   NodeId id,
                   const std::unordered_map<NodeId, std::size_t>& lastUse,
                   std::size_t i) {
    if (isImported(id)) {
      return; // imported values have no operands in this stage
    }
    std::vector<NodeId> kids;
    program_.arena().children(id, kids);
    for (const NodeId k : kids) {
      const auto it = lastUse.find(k);
      if (it != lastUse.end() && it->second == i) {
        alloc.release(k);
      }
    }
  }

  const Program& program_;
  const Analysis& analysis_;
  const ChannelMap& channels_;
  std::vector<std::int32_t>& operands_;
  std::unordered_map<NodeId, std::int32_t> boundary_;
  std::unordered_map<NodeId, std::int32_t> slotOfLiveOut_;
};

} // namespace

const char* name(Space space) {
  switch (space) {
  case Space::Uniform:
    return "uniform";
  case Space::Point:
    return "point";
  }
  return "?";
}

const char* name(Opcode op) {
  switch (op) {
  case Opcode::Const:
    return "const";
  case Opcode::LoadInput:
    return "loadin";
  case Opcode::LoadPersistent:
    return "loadp";
  case Opcode::Pw:
    return "pw";
  case Opcode::GridParam:
    return "gridparam";
  case Opcode::GridLoad:
    return "gridload";
  case Opcode::Contract:
    return "contract";
  case Opcode::Convert:
    return "convert";
  }
  return "?";
}

std::uint64_t LowerOptions::fingerprint() const {
  std::string s = "lo1;";
  std::vector<std::string> names = invariantInputs;
  std::sort(names.begin(), names.end());
  for (const std::string& n : names) {
    s += n;
    s += ',';
  }
  s += ';';
  std::vector<GridId> grids = invariantGrids;
  std::sort(grids.begin(), grids.end());
  for (const GridId g : grids) {
    s += std::to_string(g);
    s += ',';
  }
  s += ';';
  s += std::to_string(hoistThreshold);

  std::uint64_t h = FnvOffset;
  for (const char c : s) {
    h ^= static_cast<std::uint8_t>(c);
    h *= FnvPrime;
  }
  return h;
}

std::int32_t LoweredProgram::peakSlotCount() const {
  return std::max(precompute_.slotCount, run_.slotCount);
}

bool LoweredProgram::mixesTypes() const {
  for (const StageCode* stage : {&precompute_, &run_}) {
    for (const Instruction& inst : stage->code) {
      if (inst.type != computeType_) {
        return true;
      }
    }
  }
  return false;
}

std::string LoweredProgram::summary() const {
  const std::size_t hoisted = persistentSlotCount_ >= stateSlotCount_
                                  ? static_cast<std::size_t>(persistentSlotCount_ - stateSlotCount_)
                                  : 0;
  return "expr: " + std::to_string(run_.code.size()) + " ops/point, " +
         std::to_string(peakSlotCount()) + " live slots, " + std::to_string(hoisted) +
         " values hoisted out of the call (" + std::to_string(precompute_.code.size()) +
         " precompute ops), " + std::to_string(stateSlotCount_) + " state slots";
}

std::string LoweredProgram::dump() const {
  std::string out;
  const auto stage = [&](const char* title, const StageCode& s) {
    out += title;
    out += " (slots=" + std::to_string(s.slotCount) + ")\n";
    for (const Instruction& i : s.code) {
      const char* types[] = {":f32", ":f64", ":i64", ":bool"};
      out += "  t" + std::to_string(i.dst) + types[static_cast<std::size_t>(i.type)] + " = " +
             name(i.op);
      if (i.op == Opcode::Pw) {
        out += std::string(" ") + name(i.fn);
      }
      if (i.op == Opcode::Const) {
        std::array<char, 32> buf{};
        std::snprintf(buf.data(), buf.size(), "%.17g", i.value);
        out += " ";
        out += buf.data();
      }
      if (i.op == Opcode::LoadInput || i.op == Opcode::LoadPersistent) {
        out += " #" + std::to_string(i.imm);
      }
      if (i.op == Opcode::GridParam) {
        out += " g" + std::to_string(i.grid) + "." + std::to_string(i.comp) + "[" +
               std::to_string(i.imm) + "]";
      }
      if (i.op == Opcode::GridLoad) {
        out += " g" + std::to_string(i.grid) + ":" + std::to_string(i.comp);
      }
      if (i.strict) {
        out += " strict";
      }
      if (i.op == Opcode::Contract) {
        out += " m" + std::to_string(i.matrix) + " b" + std::to_string(i.block);
      }
      for (std::int32_t k = 0; k < i.operandCount; ++k) {
        out += " t" + std::to_string(operands_[i.operandBegin + k]);
      }
      out += "  ; node " + std::to_string(i.node) + " " + name(i.space) + "\n";
    }
    for (const Store& st : s.outputs) {
      out += "  out[" + std::to_string(st.target) + "] <- t" + std::to_string(st.source) + "\n";
    }
    for (const Store& st : s.persistent) {
      out += "  p[" + std::to_string(st.target) + "] <- t" + std::to_string(st.source) + "\n";
    }
    for (const CheckStore& st : s.checks) {
      out += "  check g" + std::to_string(st.grid) +
             (st.kind == CheckKind::Window ? " window" : " outside") + " <- t" +
             std::to_string(st.source) + "\n";
    }
  };
  stage("precompute", precompute_);
  stage("run", run_);
  return out;
}

bool LoweredProgram::usesGrids() const {
  return std::any_of(
      gridUses_.begin(), gridUses_.end(), [](const GridUse& use) { return use.dimensions > 0; });
}

bool LoweredProgram::hasStrict() const {
  for (const StageCode* stage : {&precompute_, &run_}) {
    for (const Instruction& inst : stage->code) {
      if (inst.strict) {
        return true;
      }
    }
  }
  return false;
}

LoweredProgram lower(const Program& untyped, const LowerOptions& options) {
  // every value with its type, every literal once per type it is used in, and the lookups written
  // out
  const ExpandedProgram expanded = expandLookups(typed(untyped));
  const Program& program = expanded.program;
  const Arena& arena = program.arena();
  const ChannelMap channels(program);
  const Analysis an = analyse(program, options, channels, expanded.checks);
  // The type the hoisted values are kept in: a value in binary64 or a whole number is hoisted only
  // where the data moves in binary64, so that keeping it loses nothing -- e.g. the weights and the
  // offsets of a lookup, which are so also in a program of binary32.
  const ComputeType data = dataType(untyped);

  // Backward sweep from the stores. A value that is invariant and expensive
  // enough becomes a stage boundary and is not descended into; everything else
  // is needed in the run stage and pulls its operands along. Invariant but cheap
  // values therefore land in the run stage and are simply recomputed, which is
  // what we want — materialising them costs a persistent buffer and a streaming
  // read to save an add.
  std::vector<char> neededInRun(arena.size(), 0);
  std::vector<char> hoisted(arena.size(), 0);
  std::vector<NodeId> stack;
  std::vector<NodeId> kids;

  const auto require = [&](NodeId id) {
    if (id == NoNode) {
      fail("a root is unset");
    }
    // Uniform-space values are deliberately excluded. They are scalars, so
    // materialising one would spend numPoints * sizeof(ComputeType) of
    // persistent memory to cache a number; the right treatment for a
    // constant-only subtree is constant folding at the frontend, not hoisting.
    const ComputeType type = resolve(arena[id].type, program.computeType());
    const bool keepable =
        data == ComputeType::F64 || type == ComputeType::F32 || type == ComputeType::Bool;
    if (an.invariant[id] != 0 && an.space[id] == Space::Point &&
        an.cost[id] >= options.hoistThreshold && keepable) {
      hoisted[id] = 1;
      return;
    }
    if (neededInRun[id] == 0) {
      neededInRun[id] = 1;
      stack.push_back(id);
    }
  };

  for (const NodeId r : program.roots()) {
    require(r);
  }
  for (const StateSpec& s : program.state()) {
    require(s.root);
  }
  for (const GridCheck& check : expanded.checks) {
    require(check.root);
  }
  while (!stack.empty()) {
    const NodeId id = stack.back();
    stack.pop_back();
    arena.children(id, kids);
    for (const NodeId k : kids) {
      require(k);
    }
  }

  // Precompute closure: everything reachable from a hoisted value. All invariant
  // by construction, so the closure never pulls a run-stage value in.
  std::vector<char> inPrecompute(arena.size(), 0);
  for (NodeId id = 0; id < static_cast<NodeId>(arena.size()); ++id) {
    if (hoisted[id] != 0) {
      stack.push_back(id);
    }
  }
  while (!stack.empty()) {
    const NodeId id = stack.back();
    stack.pop_back();
    if (inPrecompute[id] != 0) {
      continue;
    }
    inPrecompute[id] = 1;
    arena.children(id, kids);
    for (const NodeId k : kids) {
      if (inPrecompute[k] == 0) {
        stack.push_back(k);
      }
    }
  }

  LoweredProgram out;
  out.stateSlotCount_ = static_cast<std::int32_t>(program.state().size());
  out.computeType_ = program.computeType();
  out.dataType_ = data;
  out.gridUses_ = expanded.uses;

  // Persistent slot layout: declared state first, so Binding can initialise
  // [0, stateSlotCount) straight from StateSpec::initial, then hoisted values in
  // ascending NodeId order so the assignment is deterministic.
  std::unordered_map<NodeId, std::int32_t> persistentSlot;
  std::int32_t nextPersistent = out.stateSlotCount_;
  std::vector<NodeId> hoistedOrder;
  for (NodeId id = 0; id < static_cast<NodeId>(arena.size()); ++id) {
    if (hoisted[id] != 0) {
      persistentSlot.emplace(id, nextPersistent++);
      hoistedOrder.push_back(id);
    }
  }
  out.persistentSlotCount_ = nextPersistent;

  Emitter emitter(program, an, channels, out.operands_);

  if (!hoistedOrder.empty()) {
    std::vector<NodeId> values;
    for (NodeId id = 0; id < static_cast<NodeId>(arena.size()); ++id) {
      if (inPrecompute[id] != 0) {
        values.push_back(id);
      }
    }
    emitter.setBoundary({});
    emitter.emit(values, hoistedOrder, out.precompute_);
    for (const NodeId id : hoistedOrder) {
      out.precompute_.persistent.push_back(
          Store{persistentSlot.at(id), emitter.liveOutSlot(id), an.space[id]});
    }
  }

  {
    std::vector<NodeId> values;
    for (NodeId id = 0; id < static_cast<NodeId>(arena.size()); ++id) {
      if (neededInRun[id] != 0 || hoisted[id] != 0) {
        values.push_back(id);
      }
    }
    std::vector<NodeId> liveOut;
    for (const NodeId r : program.roots()) {
      liveOut.push_back(r);
    }
    for (const StateSpec& s : program.state()) {
      liveOut.push_back(s.root);
    }
    for (const GridCheck& check : expanded.checks) {
      liveOut.push_back(check.root);
    }
    emitter.setBoundary(persistentSlot);
    emitter.emit(values, liveOut, out.run_);

    for (std::size_t i = 0; i < program.roots().size(); ++i) {
      const NodeId r = program.roots()[i];
      out.run_.outputs.push_back(
          Store{static_cast<std::int32_t>(i), emitter.liveOutSlot(r), an.space[r]});
    }
    for (std::size_t i = 0; i < program.state().size(); ++i) {
      const NodeId r = program.state()[i].root;
      out.run_.persistent.push_back(
          Store{static_cast<std::int32_t>(i), emitter.liveOutSlot(r), an.space[r]});
    }
    for (const GridCheck& check : expanded.checks) {
      out.run_.checks.push_back(
          CheckStore{check.grid, check.kind, emitter.liveOutSlot(check.root)});
    }
  }

  return out;
}

} // namespace seissol::expr
