// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/Program.h"

#include "Expr/Ir.h"
#include "utils/logger.h"

#include <cstdint>
#include <cstring>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace seissol::expr {

namespace {

constexpr std::uint64_t FnvOffset = 0xCBF29CE484222325ULL;
constexpr std::uint64_t FnvPrime = 0x100000001B3ULL;

std::uint64_t fnv1a(const std::string& s) {
  std::uint64_t h = FnvOffset;
  for (const char c : s) {
    h ^= static_cast<std::uint8_t>(c);
    h *= FnvPrime;
  }
  return h;
}

void appendHex(std::string& out, std::uint64_t v) {
  static constexpr char Digits[] = "0123456789abcdef";
  for (int shift = 60; shift >= 0; shift -= 4) {
    out.push_back(Digits[(v >> static_cast<unsigned>(shift)) & 0xFU]);
  }
}

void appendInt(std::string& out, std::int64_t v) { out += std::to_string(v); }

// Const is serialised by bit pattern for the same reason interning keys on it:
// a decimal spelling makes the formatting routine part of the program identity,
// and NaN payloads and signed zero survive neither round trip reliably.
std::uint64_t bits(double v) {
  std::uint64_t u = 0;
  std::memcpy(&u, &v, sizeof(u));
  return u;
}

// Resolves a Field's interned channel id to the index the emitter sees, which is
// either the input slot the gather writes or the persistent state slot.
class ChannelResolver {
  public:
  explicit ChannelResolver(const Program& program) {
    const Arena& arena = program.arena();
    for (std::size_t i = 0; i < program.inputs().size(); ++i) {
      const int ch = arena.findChannel(program.inputs()[i].name);
      if (ch >= 0) {
        inputIndex_.emplace(ch, static_cast<int>(i));
      }
    }
    for (std::size_t i = 0; i < program.state().size(); ++i) {
      const int ch = arena.findChannel(program.state()[i].name);
      if (ch >= 0) {
        stateIndex_.emplace(ch, static_cast<int>(i));
      }
    }
  }

  // Returns 'i' or 's' in `tag` and the index; false when the channel is bound
  // to neither, which validate() rejects.
  bool resolve(int channel, char& tag, int& index) const {
    const auto in = inputIndex_.find(channel);
    if (in != inputIndex_.end()) {
      tag = 'i';
      index = in->second;
      return true;
    }
    const auto st = stateIndex_.find(channel);
    if (st != stateIndex_.end()) {
      tag = 's';
      index = st->second;
      return true;
    }
    return false;
  }

  private:
  std::unordered_map<int, int> inputIndex_;
  std::unordered_map<int, int> stateIndex_;
};

// Deterministic depth-first post-order. The order depends only on the DAG shape
// and the root order, never on NodeId values, so two frontends that build the
// same expression in different construction orders agree.
class CanonicalWalk {
  public:
  CanonicalWalk(const Program& program, const ChannelResolver& resolver)
      : program_(program), resolver_(resolver), index_(program.arena().size(), -1) {}

  int visit(NodeId root) {
    // Explicit stack: machine-generated DAGs are deep enough that recursion here
    // is a stack-overflow waiting to happen (the sderiv walks use std::function
    // recursion and inherit that risk).
    struct Frame {
      NodeId node;
      std::size_t next;
    };
    std::vector<Frame> stack;
    std::vector<NodeId> kids;
    stack.push_back({root, 0});
    while (!stack.empty()) {
      Frame& top = stack.back();
      if (index_[top.node] >= 0) {
        stack.pop_back();
        continue;
      }
      program_.arena().children(top.node, kids);
      if (top.next < kids.size()) {
        const NodeId child = kids[top.next];
        ++top.next;
        if (index_[child] < 0) {
          stack.push_back({child, 0});
        }
        continue;
      }
      emit(top.node);
      stack.pop_back();
    }
    return index_[root];
  }

  [[nodiscard]] const std::string& text() const { return text_; }

  private:
  void emit(NodeId id) {
    const Arena& arena = program_.arena();
    const Node& n = arena[id];
    const int slot = counter_++;
    index_[id] = slot;

    appendInt(text_, slot);
    text_ += '=';
    text_ += name(n.kind);
    switch (n.kind) {
    case Kind::Const:
      text_ += ' ';
      appendHex(text_, bits(n.value));
      break;
    case Kind::Field: {
      char tag = '?';
      int idx = -1;
      if (!resolver_.resolve(n.ch, tag, idx)) {
        logError() << "expr: channel" << arena.channelName(n.ch)
                   << "is neither an input nor a state; run validate() first";
      }
      text_ += ' ';
      text_ += tag;
      appendInt(text_, idx);
      break;
    }
    case Kind::PW:
      text_ += ' ';
      text_ += name(n.fn);
      appendOperands(id);
      break;
    case Kind::Lookup:
      text_ += ' ';
      appendInt(text_, n.grid);
      text_ += ':';
      appendInt(text_, n.comp);
      appendOperands(id);
      break;
    case Kind::Contract:
      text_ += ' ';
      appendInt(text_, n.matrix);
      text_ += ':';
      appendInt(text_, n.block);
      break;
    case Kind::Dx:
      text_ += ' ';
      appendInt(text_, n.axis);
      appendOperands(id);
      break;
    case Kind::Fold:
      text_ += ' ';
      text_ += name(n.red);
      appendOperands(id);
      break;
    case Kind::Cumint:
    case Kind::Sample:
      appendOperands(id);
      break;
    }
    text_ += '\n';
  }

  void appendOperands(NodeId id) {
    std::vector<NodeId> kids;
    program_.arena().children(id, kids);
    for (const NodeId k : kids) {
      text_ += ' ';
      appendInt(text_, index_[k]);
    }
  }

  const Program& program_;
  const ChannelResolver& resolver_;
  std::vector<int> index_;
  std::string text_;
  int counter_{0};
};

} // namespace

std::string Program::canonicalForm() const {
  const ChannelResolver resolver(*this);
  CanonicalWalk walk(*this, resolver);

  std::string header = "expr1 ct";
  header += (computeType_ == ComputeType::F32 ? '4' : '8');
  header += '\n';
  for (std::size_t i = 0; i < grids_.size(); ++i) {
    header += "grid ";
    appendInt(header, static_cast<std::int64_t>(i));
    header += ' ';
    header += grids_[i].canonicalKey();
    header += '\n';
  }
  // Matrices and blocks by index and form only: their names reach the binding, not the emitter.
  for (std::size_t i = 0; i < matrices_.size(); ++i) {
    header += "matrix ";
    appendInt(header, static_cast<std::int64_t>(i));
    header += ' ';
    appendInt(header, static_cast<std::int64_t>(matrices_[i].shape.rows));
    header += 'x';
    appendInt(header, static_cast<std::int64_t>(matrices_[i].shape.cols));
    header += '/';
    appendInt(header, static_cast<std::int64_t>(matrices_[i].shape.leadingDimension));
    header += '\n';
  }
  for (std::size_t i = 0; i < blocks_.size(); ++i) {
    header += "block ";
    appendInt(header, static_cast<std::int64_t>(i));
    header += ' ';
    appendInt(header, static_cast<std::int64_t>(blocks_[i].length));
    header += '\n';
  }

  // Roots first, then state roots, both in declaration order: the two lists are
  // separate targets and swapping one for the other is a different kernel.
  std::string tail;
  for (const NodeId root : roots_) {
    tail += "out ";
    appendInt(tail, walk.visit(root));
    tail += '\n';
  }
  for (const StateSpec& s : state_) {
    tail += "state ";
    appendInt(tail, walk.visit(s.root));
    tail += ' ';
    appendHex(tail, bits(s.initial));
    tail += '\n';
  }

  return header + walk.text() + tail;
}

std::uint64_t Program::fingerprint() const { return fnv1a(canonicalForm()); }

void Program::addOutput(const std::string& name, reader::scripting::DataType type, NodeId root) {
  outputs_.push_back(VarSpec{name, type});
  roots_.push_back(root);
}

void Program::addInput(const std::string& name, reader::scripting::DataType type) {
  inputs_.push_back(VarSpec{name, type});
}

void Program::addState(const std::string& name, double initial, NodeId root) {
  state_.push_back(StateSpec{name, initial, root});
}

MatrixId Program::internMatrix(const std::string& name, const MatrixShape& shape) {
  if (shape.rows == 0 || shape.cols == 0 || shape.leadingDimension < shape.rows) {
    throw std::invalid_argument("expr: matrix '" + name + "' has the invalid form " +
                                std::to_string(shape.rows) + "x" + std::to_string(shape.cols) +
                                " with leading dimension " +
                                std::to_string(shape.leadingDimension));
  }
  for (std::size_t i = 0; i < matrices_.size(); ++i) {
    if (matrices_[i].name == name) {
      if (matrices_[i].shape != shape) {
        throw std::invalid_argument("expr: matrix '" + name + "' is declared with two forms");
      }
      return static_cast<MatrixId>(i);
    }
  }
  if (!matrices_.empty() && matrices_.front().shape.rows != shape.rows) {
    throw std::invalid_argument("expr: matrix '" + name + "' has " + std::to_string(shape.rows) +
                                " rows, but the point set has " +
                                std::to_string(matrices_.front().shape.rows) + " points per cell");
  }
  matrices_.push_back(MatrixSpec{name, shape});
  return static_cast<MatrixId>(matrices_.size() - 1);
}

BlockId Program::internBlock(const std::string& name, std::size_t length) {
  if (length == 0) {
    throw std::invalid_argument("expr: block '" + name + "' has no coefficients");
  }
  for (std::size_t i = 0; i < blocks_.size(); ++i) {
    if (blocks_[i].name == name) {
      if (blocks_[i].length != length) {
        throw std::invalid_argument("expr: block '" + name + "' is declared with two lengths");
      }
      return static_cast<BlockId>(i);
    }
  }
  blocks_.push_back(BlockSpec{name, length});
  return static_cast<BlockId>(blocks_.size() - 1);
}

GridId Program::internGrid(const reader::datafield::GridDesc& desc) {
  const std::string key = desc.canonicalKey() + '\x1f' + desc.path + '\x1f' + desc.variable;
  for (std::size_t i = 0; i < grids_.size(); ++i) {
    const std::string other =
        grids_[i].canonicalKey() + '\x1f' + grids_[i].path + '\x1f' + grids_[i].variable;
    if (other == key) {
      return static_cast<GridId>(i);
    }
  }
  grids_.push_back(desc);
  return static_cast<GridId>(grids_.size() - 1);
}

std::set<std::string> channelsRead(const Program& program, const std::vector<NodeId>& roots) {
  std::set<std::string> read;
  std::vector<bool> seen(program.arena().size(), false);
  std::vector<NodeId> pending = roots;
  std::vector<NodeId> children;
  while (!pending.empty()) {
    const auto node = pending.back();
    pending.pop_back();
    if (node == NoNode || seen[static_cast<std::size_t>(node)]) {
      continue;
    }
    seen[static_cast<std::size_t>(node)] = true;
    if (program.arena()[node].kind == Kind::Field) {
      read.insert(program.arena().channelName(program.arena()[node].ch));
    }
    program.arena().children(node, children);
    pending.insert(pending.end(), children.begin(), children.end());
  }
  return read;
}

} // namespace seissol::expr
