// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/Ir.h"

#include "utils/logger.h"

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

namespace seissol::expr {

namespace {

constexpr std::uint64_t HashSeed = 0x9E3779B97F4A7C15ULL;

std::size_t mix(std::size_t h, std::size_t v) { return h ^ (v + HashSeed + (h << 6U) + (h >> 2U)); }

// Bit pattern of a double, not its value. This is the whole point of delta 4:
// `==` is not reflexive for NaN, so a NaN key inserted into the pool can never
// be found again and every konst(NaN) allocates a fresh node. It is also not
// discriminating enough for signed zero — +0.0 == -0.0 is true, so the two
// collapse to one node and 1.0/x silently changes sign for whichever of the two
// lost the race. Bit patterns fix both directions at once.
std::uint64_t bits(double v) {
  std::uint64_t u = 0;
  std::memcpy(&u, &v, sizeof(u));
  return u;
}

} // namespace

const char* name(Kind kind) {
  switch (kind) {
  case Kind::Const:
    return "const";
  case Kind::Field:
    return "field";
  case Kind::PW:
    return "pw";
  case Kind::Lookup:
    return "lookup";
  case Kind::Dx:
    return "dx";
  case Kind::Cumint:
    return "cumint";
  case Kind::Fold:
    return "fold";
  case Kind::Sample:
    return "sample";
  }
  return "?";
}

bool NodeEq::operator()(const Node& x, const Node& y) const {
  if (x.kind != y.kind) {
    return false;
  }
  switch (x.kind) {
  case Kind::Const:
    return bits(x.value) == bits(y.value);
  case Kind::Field:
    return x.ch == y.ch;
  case Kind::PW: {
    if (x.fn != y.fn) {
      return false;
    }
    const int n = arity(x.fn);
    return x.a == y.a && (n < 2 || x.b == y.b) && (n < 3 || x.c == y.c);
  }
  case Kind::Lookup: {
    if (x.grid != y.grid || x.comp != y.comp || x.argCount != y.argCount) {
      return false;
    }
    const NodeId* xa = arena->args(x);
    const NodeId* ya = arena->args(y);
    for (std::int32_t i = 0; i < x.argCount; ++i) {
      if (xa[i] != ya[i]) {
        return false;
      }
    }
    return true;
  }
  case Kind::Dx:
    return x.axis == y.axis && x.a == y.a;
  case Kind::Cumint:
  case Kind::Sample:
    return x.a == y.a;
  case Kind::Fold:
    return x.red == y.red && x.a == y.a;
  }
  return false;
}

std::size_t NodeHash::operator()(const Node& n) const {
  std::size_t h = static_cast<std::size_t>(n.kind) + 1;
  switch (n.kind) {
  case Kind::Const:
    h = mix(h, static_cast<std::size_t>(bits(n.value)));
    break;
  case Kind::Field:
    h = mix(h, static_cast<std::size_t>(n.ch));
    break;
  case Kind::PW: {
    h = mix(h, static_cast<std::size_t>(n.fn));
    const int a = arity(n.fn);
    h = mix(h, static_cast<std::size_t>(n.a));
    if (a >= 2) {
      h = mix(h, static_cast<std::size_t>(n.b));
    }
    if (a >= 3) {
      h = mix(h, static_cast<std::size_t>(n.c));
    }
    break;
  }
  case Kind::Lookup: {
    h = mix(h, static_cast<std::size_t>(n.grid));
    h = mix(h, static_cast<std::size_t>(n.comp));
    h = mix(h, static_cast<std::size_t>(n.argCount));
    const NodeId* a = arena->args(n);
    for (std::int32_t i = 0; i < n.argCount; ++i) {
      h = mix(h, static_cast<std::size_t>(a[i]));
    }
    break;
  }
  case Kind::Dx:
    h = mix(h, static_cast<std::size_t>(n.axis));
    h = mix(h, static_cast<std::size_t>(n.a));
    break;
  case Kind::Cumint:
  case Kind::Sample:
    h = mix(h, static_cast<std::size_t>(n.a));
    break;
  case Kind::Fold:
    h = mix(h, static_cast<std::size_t>(n.red));
    h = mix(h, static_cast<std::size_t>(n.a));
    break;
  }
  return h;
}

Arena::Arena() : pool_(0, NodeHash{this}, NodeEq{this}) {}

Arena::Arena(const Arena& other)
    : nodes_(other.nodes_), args_(other.args_), pool_(0, NodeHash{this}, NodeEq{this}),
      channelNames_(other.channelNames_), channelIds_(other.channelIds_) {
  rebuildPool();
}

Arena& Arena::operator=(const Arena& other) {
  if (this != &other) {
    nodes_ = other.nodes_;
    args_ = other.args_;
    channelNames_ = other.channelNames_;
    channelIds_ = other.channelIds_;
    rebuildPool();
  }
  return *this;
}

Arena::Arena(Arena&& other) noexcept
    : nodes_(std::move(other.nodes_)), args_(std::move(other.args_)),
      pool_(0, NodeHash{this}, NodeEq{this}), channelNames_(std::move(other.channelNames_)),
      channelIds_(std::move(other.channelIds_)) {
  rebuildPool();
  other.pool_.clear();
}

Arena& Arena::operator=(Arena&& other) noexcept {
  if (this != &other) {
    nodes_ = std::move(other.nodes_);
    args_ = std::move(other.args_);
    channelNames_ = std::move(other.channelNames_);
    channelIds_ = std::move(other.channelIds_);
    rebuildPool();
    other.pool_.clear();
  }
  return *this;
}

// pool_ cannot simply be carried across: its hasher and comparator hold a
// back-pointer to the owning Arena, and unordered_map offers no way to re-aim
// them after construction. Rebuilding from nodes_ is O(n) and happens once per
// Program handover, which is not a hot path.
void Arena::rebuildPool() {
  pool_ = std::unordered_map<Node, NodeId, NodeHash, NodeEq>(
      nodes_.size(), NodeHash{this}, NodeEq{this});
  for (std::size_t i = 0; i < nodes_.size(); ++i) {
    pool_.emplace(nodes_[i], static_cast<NodeId>(i));
  }
}

int Arena::channel(const std::string& name) {
  const auto it = channelIds_.find(name);
  if (it != channelIds_.end()) {
    return it->second;
  }
  const int id = static_cast<int>(channelNames_.size());
  channelIds_.emplace(name, id);
  channelNames_.push_back(name);
  return id;
}

int Arena::findChannel(const std::string& name) const {
  const auto it = channelIds_.find(name);
  return it == channelIds_.end() ? -1 : it->second;
}

NodeId Arena::intern(const Node& n) {
  const auto it = pool_.find(n);
  if (it != pool_.end()) {
    return it->second;
  }
  const auto id = static_cast<NodeId>(nodes_.size());
  nodes_.push_back(n);
  pool_.emplace(n, id);
  return id;
}

NodeId Arena::konst(double v) {
  Node n;
  n.kind = Kind::Const;
  n.value = v;
  return intern(n);
}

NodeId Arena::field(int ch) {
  Node n;
  n.kind = Kind::Field;
  n.ch = ch;
  return intern(n);
}

NodeId Arena::pw(Fn f, NodeId x) {
  if (arity(f) != 1) {
    logError() << "expr: arity mismatch for" << name(f) << "— expected" << arity(f)
               << "arguments, got 1";
  }
  Node n;
  n.kind = Kind::PW;
  n.fn = f;
  n.a = x;
  return intern(n);
}

NodeId Arena::pw(Fn f, NodeId x, NodeId y) {
  if (arity(f) != 2) {
    logError() << "expr: arity mismatch for" << name(f) << "— expected" << arity(f)
               << "arguments, got 2";
  }
  Node n;
  n.kind = Kind::PW;
  n.fn = f;
  n.a = x;
  n.b = y;
  return intern(n);
}

NodeId Arena::pw(Fn f, NodeId x, NodeId y, NodeId z) {
  if (arity(f) != 3) {
    logError() << "expr: arity mismatch for" << name(f) << "— expected" << arity(f)
               << "arguments, got 3";
  }
  Node n;
  n.kind = Kind::PW;
  n.fn = f;
  n.a = x;
  n.b = y;
  n.c = z;
  return intern(n);
}

// The coordinate span has to be visible to the hasher *before* the pool lookup,
// so the arguments are appended first and rolled back on a hit. Without that,
// two structurally equal Lookups would hash differently and interning would
// silently stop working for exactly the node kind where it matters most.
NodeId Arena::lookup(GridId grid, std::int32_t component, const std::vector<NodeId>& coords) {
  if (coords.empty() || coords.size() > static_cast<std::size_t>(MaxLookupDimension)) {
    logError() << "expr: lookup needs 1 .." << MaxLookupDimension << "coordinates, got"
               << coords.size();
  }
  const std::size_t mark = args_.size();
  args_.insert(args_.end(), coords.begin(), coords.end());

  Node n;
  n.kind = Kind::Lookup;
  n.grid = grid;
  n.comp = component;
  n.argBegin = static_cast<std::int32_t>(mark);
  n.argCount = static_cast<std::int32_t>(coords.size());

  const auto it = pool_.find(n);
  if (it != pool_.end()) {
    args_.resize(mark);
    return it->second;
  }
  const auto id = static_cast<NodeId>(nodes_.size());
  nodes_.push_back(n);
  pool_.emplace(n, id);
  return id;
}

NodeId Arena::dx(int axis, NodeId x) {
  Node n;
  n.kind = Kind::Dx;
  n.axis = axis;
  n.a = x;
  return intern(n);
}

NodeId Arena::cumint(NodeId x) {
  Node n;
  n.kind = Kind::Cumint;
  n.a = x;
  return intern(n);
}

NodeId Arena::fold(Red r, NodeId x) {
  Node n;
  n.kind = Kind::Fold;
  n.red = r;
  n.a = x;
  return intern(n);
}

NodeId Arena::sample(NodeId x) {
  Node n;
  n.kind = Kind::Sample;
  n.a = x;
  return intern(n);
}

void Arena::children(NodeId id, std::vector<NodeId>& out) const {
  out.clear();
  const Node& n = nodes_[id];
  switch (n.kind) {
  case Kind::Const:
  case Kind::Field:
    return;
  case Kind::PW: {
    const int a = arity(n.fn);
    out.push_back(n.a);
    if (a >= 2) {
      out.push_back(n.b);
    }
    if (a >= 3) {
      out.push_back(n.c);
    }
    return;
  }
  case Kind::Lookup: {
    const NodeId* a = args(n);
    out.insert(out.end(), a, a + n.argCount);
    return;
  }
  case Kind::Dx:
  case Kind::Cumint:
  case Kind::Fold:
  case Kind::Sample:
    out.push_back(n.a);
    return;
  }
}

std::vector<NodeId> Arena::children(NodeId id) const {
  std::vector<NodeId> out;
  children(id, out);
  return out;
}

const char* name(Fn f) {
  switch (f) {
  case Fn::Neg:
    return "neg";
  case Fn::Sqrt:
    return "sqrt";
  case Fn::Abs:
    return "abs";
  case Fn::Exp:
    return "exp";
  case Fn::Log:
    return "log";
  case Fn::Log2:
    return "log2";
  case Fn::Log10:
    return "log10";
  case Fn::Sign:
    return "sign";
  case Fn::Floor:
    return "floor";
  case Fn::Ceil:
    return "ceil";
  case Fn::Round:
    return "round";
  case Fn::Rcp:
    return "rcp";
  case Fn::Sin:
    return "sin";
  case Fn::Cos:
    return "cos";
  case Fn::Tan:
    return "tan";
  case Fn::Asin:
    return "asin";
  case Fn::Acos:
    return "acos";
  case Fn::Atan:
    return "atan";
  case Fn::Sinh:
    return "sinh";
  case Fn::Cosh:
    return "cosh";
  case Fn::Tanh:
    return "tanh";
  case Fn::Asinh:
    return "asinh";
  case Fn::Acosh:
    return "acosh";
  case Fn::Atanh:
    return "atanh";
  case Fn::Erf:
    return "erf";
  case Fn::Add:
    return "add";
  case Fn::Sub:
    return "sub";
  case Fn::Mul:
    return "mul";
  case Fn::Div:
    return "div";
  case Fn::Pow:
    return "pow";
  case Fn::Min:
    return "min";
  case Fn::Max:
    return "max";
  case Fn::Atan2:
    return "atan2";
  case Fn::Mod:
    return "mod";
  case Fn::Lt:
    return "lt";
  case Fn::Le:
    return "le";
  case Fn::Eq:
    return "eq";
  case Fn::And:
    return "and";
  case Fn::Or:
    return "or";
  case Fn::Select:
    return "select";
  }
  return "?";
}

const char* name(Red r) {
  switch (r) {
  case Red::Max:
    return "max";
  case Red::Min:
    return "min";
  case Red::Int:
    return "int";
  case Red::Sum:
    return "sum";
  case Red::Mean:
    return "mean";
  case Red::Rms:
    return "rms";
  case Red::ArgMax:
    return "argmax";
  case Red::ArgMin:
    return "argmin";
  case Red::Last:
    return "last";
  }
  return "?";
}

namespace {

// Driven off name(Fn) so the two tables cannot drift apart. The Fn enum is dense
// and Select is its last member, which is the only ordering fact relied on here.
template <typename Enum, Enum Last>
bool fromNameImpl(const std::string& s, Enum& out) {
  for (int i = 0; i <= static_cast<int>(Last); ++i) {
    const auto e = static_cast<Enum>(i);
    if (s == name(e)) {
      out = e;
      return true;
    }
  }
  return false;
}

} // namespace

bool fnFromName(const std::string& s, Fn& out) { return fromNameImpl<Fn, Fn::Select>(s, out); }

bool redFromName(const std::string& s, Red& out) { return fromNameImpl<Red, Red::Last>(s, out); }

} // namespace seissol::expr
