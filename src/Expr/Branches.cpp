// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/Branches.h"

#include "Expr/Ir.h"
#include "Expr/Program.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

namespace seissol::expr {

namespace {

bool isSelect(const Node& node) { return node.kind == Kind::PW && node.fn == Fn::Select; }

/// The outputs, the states and the checks of `program`.
std::vector<NodeId> rootsOf(const Program& program) {
  std::vector<NodeId> roots = program.roots();
  for (const auto& state : program.state()) {
    roots.push_back(state.root);
  }
  for (const auto& check : program.checks()) {
    roots.push_back(check.root);
  }
  return roots;
}

/// What `roots` need of `arena`, in ascending order.
std::vector<NodeId> reachableFrom(const Arena& arena, const std::vector<NodeId>& roots) {
  std::vector<char> seen(arena.size(), 0);
  std::vector<NodeId> stack;
  std::vector<NodeId> kids;
  for (const NodeId root : roots) {
    if (root != NoNode) {
      stack.push_back(root);
    }
  }
  while (!stack.empty()) {
    const NodeId id = stack.back();
    stack.pop_back();
    if (seen[static_cast<std::size_t>(id)] != 0) {
      continue;
    }
    seen[static_cast<std::size_t>(id)] = 1;
    arena.children(id, kids);
    stack.insert(stack.end(), kids.begin(), kids.end());
  }
  std::vector<NodeId> values;
  for (std::size_t id = 0; id < seen.size(); ++id) {
    if (seen[id] != 0) {
      values.push_back(static_cast<NodeId>(id));
    }
  }
  return values;
}

} // namespace

Regions::Regions(const Arena& arena,
                 const std::vector<NodeId>& values,
                 const std::vector<NodeId>& roots,
                 const std::vector<char>* leaves,
                 bool truths)
    : truths_(truths), regions_(1), regionOf_(arena.size(), -1) {
  for (const NodeId root : roots) {
    if (root != NoNode) {
      meet(root, 0);
    }
  }
  // the users of a value before it: its uses are known when it is reached
  std::vector<NodeId> kids;
  for (auto it = values.rbegin(); it != values.rend(); ++it) {
    const NodeId value = *it;
    const std::int32_t region = regionOf_[static_cast<std::size_t>(value)];
    if (region < 0 || (leaves != nullptr && (*leaves)[static_cast<std::size_t>(value)] != 0)) {
      continue;
    }
    const Node& node = arena[value];
    arena.children(value, kids);
    for (std::size_t k = 0; k < kids.size(); ++k) {
      meet(kids[k], k > 0 && branches(arena, node) ? branch(region, node.a, k == 1) : region);
    }
  }
}

std::int32_t Regions::of(NodeId value) const {
  return value >= 0 && static_cast<std::size_t>(value) < regionOf_.size()
             ? regionOf_[static_cast<std::size_t>(value)]
             : -1;
}

std::int32_t Regions::common(std::int32_t x, std::int32_t y) const {
  while ((*this)[x].depth > (*this)[y].depth) {
    x = (*this)[x].parent;
  }
  while ((*this)[y].depth > (*this)[x].depth) {
    y = (*this)[y].parent;
  }
  while (x != y) {
    x = (*this)[x].parent;
    y = (*this)[y].parent;
  }
  return x;
}

std::int32_t Regions::useOf(const Arena& arena,
                            NodeId user,
                            std::int32_t region,
                            std::size_t operand) const {
  const Node& node = arena[user];
  if (operand == 0 || !branches(arena, node)) {
    return region;
  }
  const auto it = branches_.find({region, node.a, operand == 1});
  return it == branches_.end() ? region : it->second;
}

bool Regions::branches(const Arena& arena, const Node& node) const {
  return isSelect(node) && (truths_ || arena[node.a].type != ValueType::Bool);
}

std::int32_t Regions::branch(std::int32_t parent, NodeId condition, bool side) {
  const auto key = std::make_tuple(parent, condition, side);
  const auto it = branches_.find(key);
  if (it != branches_.end()) {
    return it->second;
  }
  const auto index = static_cast<std::int32_t>(regions_.size());
  regions_.push_back(Region{parent, condition, side, (*this)[parent].depth + 1});
  branches_.emplace(key, index);
  return index;
}

void Regions::meet(NodeId value, std::int32_t region) {
  std::int32_t& current = regionOf_[static_cast<std::size_t>(value)];
  current = current < 0 ? region : common(current, region);
}

Reach::Reach(const Program& program,
             Arena& target,
             std::function<NodeId(NodeId)> condition,
             ValueType type)
    : arena_(program.arena()), target_(target), condition_(std::move(condition)), type_(type),
      roots_(rootsOf(program)), values_(reachableFrom(arena_, roots_)),
      regions_(arena_, values_, roots_), uses_(arena_.size()), done_(arena_.size(), 0),
      within_(arena_.size(), NoNode), plain_(arena_.size(), -1) {
  std::vector<NodeId> kids;
  for (const NodeId value : values_) {
    const std::int32_t region = regions_.of(value);
    if (region < 0) {
      continue;
    }
    arena_.children(value, kids);
    for (std::size_t k = 0; k < kids.size(); ++k) {
      uses_[static_cast<std::size_t>(kids[k])].push_back(
          Use{value, regions_.useOf(arena_, value, region, k)});
    }
  }
  for (const NodeId root : roots_) {
    if (root != NoNode) {
      uses_[static_cast<std::size_t>(root)].push_back(Use{NoNode, 0});
    }
  }
  true_ = target_.konst(1.0, type_);
  false_ = target_.konst(0.0, type_);
}

NodeId Reach::choose(NodeId condition, bool side, NodeId x, NodeId y) {
  if (x == y) {
    return x;
  }
  return side ? target_.pw(Fn::Select, condition, x, y, type_)
              : target_.pw(Fn::Select, condition, y, x, type_);
}

NodeId Reach::masked(NodeId node, NodeId value, NodeId otherwise) {
  const std::int32_t region = regions_.of(node);
  if (region < 0) {
    return otherwise;
  }
  NodeId result = value;
  const NodeId where = within(node);
  if (where != NoNode) {
    result = choose(where, true, result, otherwise);
  }
  for (std::int32_t r = region; r > 0; r = regions_[r].parent) {
    result = choose(condition_(regions_[r].condition), regions_[r].side, result, otherwise);
  }
  return result;
}

bool Reach::everywhere(NodeId node) {
  if (regions_.of(node) != 0) {
    return false;
  }
  // used everywhere in its region: by a root, or in that region by a value used everywhere in it,
  // as withinFromUses() finds -- the users first, without recursion
  const auto plain = [&](NodeId value) -> signed char& {
    return plain_[static_cast<std::size_t>(value)];
  };
  std::vector<NodeId> stack{node};
  while (!stack.empty()) {
    const NodeId value = stack.back();
    if (plain(value) >= 0) {
      stack.pop_back();
      continue;
    }
    const std::int32_t home = regions_.of(value);
    bool pending = false;
    for (const Use& use : uses_[static_cast<std::size_t>(value)]) {
      if (use.user != NoNode && use.region == home && plain(use.user) < 0) {
        stack.push_back(use.user);
        pending = true;
      }
    }
    if (pending) {
      continue;
    }
    stack.pop_back();
    plain(value) = 0;
    for (const Use& use : uses_[static_cast<std::size_t>(value)]) {
      if (use.region == home && (use.user == NoNode || plain(use.user) == 1)) {
        plain(value) = 1;
        break;
      }
    }
  }
  return plain(node) == 1;
}

NodeId Reach::within(NodeId node) {
  // the users first, without recursion: a chain of uses is as long as the program is deep
  std::vector<NodeId> stack{node};
  while (!stack.empty()) {
    const NodeId value = stack.back();
    if (done_[static_cast<std::size_t>(value)] != 0) {
      stack.pop_back();
      continue;
    }
    bool pending = false;
    for (const Use& use : uses_[static_cast<std::size_t>(value)]) {
      if (use.user != NoNode && done_[static_cast<std::size_t>(use.user)] == 0) {
        stack.push_back(use.user);
        pending = true;
      }
    }
    if (pending) {
      continue;
    }
    stack.pop_back();
    within_[static_cast<std::size_t>(value)] = withinFromUses(value);
    done_[static_cast<std::size_t>(value)] = 1;
  }
  return within_[static_cast<std::size_t>(node)];
}

NodeId Reach::withinFromUses(NodeId node) {
  const std::int32_t home = regions_.of(node);
  std::vector<NodeId> terms;
  for (const Use& use : uses_[static_cast<std::size_t>(node)]) {
    // where the user is used in its region, and the branches between that region and `home`
    const NodeId user = use.user == NoNode ? NoNode : within_[static_cast<std::size_t>(use.user)];
    if (use.region == home && user == NoNode) {
      return NoNode;
    }
    NodeId term = user == NoNode ? true_ : user;
    for (std::int32_t r = use.region; r != home; r = regions_[r].parent) {
      term = choose(condition_(regions_[r].condition), regions_[r].side, term, false_);
    }
    if (std::find(terms.begin(), terms.end(), term) == terms.end()) {
      terms.push_back(term);
    }
  }
  if (terms.empty()) {
    return false_;
  }
  // either of them, the first tested first
  NodeId result = terms.back();
  for (auto it = std::next(terms.rbegin()); it != terms.rend(); ++it) {
    result = choose(*it, true, true_, result);
  }
  return result;
}

std::vector<Step> branchedOrder(const Arena& arena,
                                const std::vector<NodeId>& values,
                                const std::vector<NodeId>& roots,
                                const std::vector<char>& leaves,
                                const std::function<std::int64_t(NodeId)>& weight,
                                std::int64_t threshold) {
  std::vector<Step> steps;
  steps.reserve(values.size());
  const auto straight = [&] {
    steps.clear();
    for (const NodeId value : values) {
      steps.push_back(Step{Step::Kind::Value, value, true});
    }
    return steps;
  };
  const Regions regions(arena, values, roots, &leaves, false);
  const auto count = static_cast<std::size_t>(regions.size());
  if (count == 1) {
    return straight();
  }

  // what the values of a region weigh, with those of the regions in it
  constexpr std::int64_t Cap = std::numeric_limits<std::int64_t>::max() / 4;
  std::vector<std::int64_t> weights(count, 0);
  for (const NodeId value : values) {
    const std::int32_t region = regions.of(value);
    if (region >= 0) {
      auto& w = weights[static_cast<std::size_t>(region)];
      w = std::min(Cap, w + weight(value));
    }
  }
  for (std::size_t r = count - 1; r > 0; --r) {
    auto& w = weights[static_cast<std::size_t>(regions[static_cast<std::int32_t>(r)].parent)];
    w = std::min(Cap, w + weights[r]);
  }
  // a region is laid out in the nearest branch around it, or in region 0: a branch is created
  // after the region it lies in
  std::vector<std::int32_t> layout(count, 0);
  std::vector<std::int32_t> outer(count, 0);
  std::vector<std::vector<std::int32_t>> inner(count);
  bool any = false;
  for (std::size_t r = 1; r < count; ++r) {
    const auto parent = static_cast<std::size_t>(regions[static_cast<std::int32_t>(r)].parent);
    if (weights[r] >= threshold) {
      layout[r] = static_cast<std::int32_t>(r);
      outer[r] = layout[parent];
      inner[static_cast<std::size_t>(outer[r])].push_back(static_cast<std::int32_t>(r));
      any = true;
    } else {
      layout[r] = layout[parent];
    }
  }
  if (!any) {
    return straight();
  }
  const auto layoutOf = [&](NodeId value) {
    const std::int32_t region = regions.of(value);
    return region < 0 ? 0 : layout[static_cast<std::size_t>(region)];
  };

  // a list schedule: of the values whose operands are computed, those of the current branch by
  // NodeId; then the branch in it of the first value ready in it; the branch ends where neither is
  // left
  std::vector<std::int32_t> remaining(arena.size(), 0);
  std::vector<std::vector<NodeId>> users(arena.size());
  std::vector<NodeId> kids;
  for (const NodeId value : values) {
    if (leaves[static_cast<std::size_t>(value)] != 0) {
      continue;
    }
    arena.children(value, kids);
    remaining[static_cast<std::size_t>(value)] = static_cast<std::int32_t>(kids.size());
    for (const NodeId kid : kids) {
      users[static_cast<std::size_t>(kid)].push_back(value);
    }
  }
  std::vector<std::set<NodeId>> ready(count);
  std::vector<std::int64_t> readyWithin(count, 0);
  std::vector<char> computed(arena.size(), 0);
  const auto mark = [&](NodeId value, std::int64_t delta) {
    for (std::int32_t r = layoutOf(value);; r = outer[static_cast<std::size_t>(r)]) {
      readyWithin[static_cast<std::size_t>(r)] += delta;
      if (r == 0) {
        break;
      }
    }
  };
  const auto makeReady = [&](NodeId value) {
    ready[static_cast<std::size_t>(layoutOf(value))].insert(value);
    mark(value, 1);
  };
  for (const NodeId value : values) {
    if (remaining[static_cast<std::size_t>(value)] == 0) {
      makeReady(value);
    }
  }

  // the first value ready in `region` or in a branch in it whose condition is computed
  const std::function<NodeId(std::int32_t)> first = [&](std::int32_t region) {
    const auto r = static_cast<std::size_t>(region);
    if (readyWithin[r] == 0) {
      return NoNode;
    }
    NodeId best = ready[r].empty() ? NoNode : *ready[r].begin();
    for (const std::int32_t branch : inner[r]) {
      if (computed[static_cast<std::size_t>(regions[branch].condition)] == 0) {
        continue;
      }
      const NodeId candidate = first(branch);
      if (candidate != NoNode && (best == NoNode || candidate < best)) {
        best = candidate;
      }
    }
    return best;
  };

  // the other side of the condition of a branch, in the same branch around: after the one, the
  // other where it has values, so that a backend makes an `else` of it
  std::map<std::tuple<std::int32_t, NodeId, bool>, std::int32_t> sides;
  for (std::size_t r = 1; r < count; ++r) {
    if (layout[r] == static_cast<std::int32_t>(r)) {
      const Region& region = regions[static_cast<std::int32_t>(r)];
      sides.emplace(std::make_tuple(outer[r], region.condition, region.side),
                    static_cast<std::int32_t>(r));
    }
  }

  std::int32_t current = 0;
  std::int32_t closed = -1;
  std::size_t done = 0;
  while (done < values.size()) {
    if (closed > 0) {
      const Region& region = regions[closed];
      const auto other = sides.find(
          std::make_tuple(outer[static_cast<std::size_t>(closed)], region.condition, !region.side));
      closed = -1;
      if (other != sides.end() && readyWithin[static_cast<std::size_t>(other->second)] > 0) {
        steps.push_back(Step{Step::Kind::If, region.condition, !region.side});
        current = other->second;
        continue;
      }
    }
    auto& here = ready[static_cast<std::size_t>(current)];
    if (here.empty()) {
      // into the branch in the current one of the first value ready in it, or out of this one
      NodeId next = NoNode;
      std::int32_t into = -1;
      for (const std::int32_t branch : inner[static_cast<std::size_t>(current)]) {
        if (computed[static_cast<std::size_t>(regions[branch].condition)] == 0) {
          continue;
        }
        const NodeId candidate = first(branch);
        if (candidate != NoNode && (next == NoNode || candidate < next)) {
          next = candidate;
          into = branch;
        }
      }
      if (into >= 0) {
        steps.push_back(Step{Step::Kind::If, regions[into].condition, regions[into].side});
        current = into;
        continue;
      }
      if (current == 0) {
        // no value is left whose branches can start: the values of a branch and a condition of
        // one around it need each other, which the arena should not allow; computed everywhere,
        // they give the same values
        return straight();
      }
      steps.push_back(Step{Step::Kind::EndIf, NoNode, true});
      closed = current;
      current = outer[static_cast<std::size_t>(current)];
      continue;
    }
    // the values of the current branch first, so that a branch it holds starts once, with what it
    // needs of this one computed
    const NodeId next = *here.begin();
    steps.push_back(Step{Step::Kind::Value, next, true});
    ready[static_cast<std::size_t>(current)].erase(next);
    mark(next, -1);
    computed[static_cast<std::size_t>(next)] = 1;
    ++done;
    for (const NodeId user : users[static_cast<std::size_t>(next)]) {
      if (--remaining[static_cast<std::size_t>(user)] == 0) {
        makeReady(user);
      }
    }
  }
  while (current != 0) {
    steps.push_back(Step{Step::Kind::EndIf, NoNode, true});
    current = outer[static_cast<std::size_t>(current)];
  }
  return steps;
}

} // namespace seissol::expr
