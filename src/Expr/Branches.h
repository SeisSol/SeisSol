// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_BRANCHES_H_
#define SEISSOL_SRC_EXPR_BRANCHES_H_

// The branches of a program (part of the lowering).
//
// A select is a value of the point: the operand its condition chooses, the first where the
// condition is not zero (a NaN included). Only that operand reaches an output there, so a value
// that one branch of a select alone uses is needed only at the points that take the branch: a grid
// sampled in one branch of a dispatch over the layers of a model is needed only in that layer. This
// finds where each value is needed, and the lowering makes two things of it: branches of the code
// around the values that only some points need, which each backend runs as its hardware wants; and
// the checks of the lookups, which count only where their sample is used.
//
// REGIONS. The points of a region R that reach a select of the condition c (the select lies in R)
// use its first operand in the region (R, c, true), where c holds, its second in (R, c, false), and
// its condition in R. A value lies in the innermost region that holds all its uses (Regions); an
// output, a state or a check of the program uses its value in region 0, at every point. The
// regions are those of conditions, not of selects: the selects a frontend makes of one `if`, one
// per variable it assigns, share their branches, and the checks of a sample share them with it.
//
// WHERE A VALUE IS USED (Reach). The region of a value holds every point where it is used, and may
// hold more: a value that branches of two conditions use lies where both lie, but is used only
// where either is taken. A check of a lookup -- a sample beyond the volume of a grid that fails
// there, or at a time its window does not hold -- counts exactly where its sample is used: where
// the selects on some way from it to an output, a state or a check of the program all take the
// branch it lies in. Its truth value is written as selects over the conditions of the regions
// around the sample (Reach::masked()), so that it lies in the region of the sample, and means the
// same on every backend, whether it runs a branch or computes it everywhere.
//
// THE ORDER (branchedOrder()). A region becomes a branch of the code where its values weigh enough
// to be worth one; the values of a lighter region are computed with those around it, at every
// point, as all values were before there were branches. So are those of the selects whose condition
// is a truth value (ValueType::Bool): those of a written-out lookup, between the cases of a
// coordinate at the ends of an axis, and those of where a sample is used (Reach). Nearly every
// point takes the same side of them, so a branch would cost its test and the order of the values
// around it, and save nothing. In a region, the values that can be computed come first, then the
// branches in it, so that a branch starts once, with what it needs of the values around it
// computed, and holds its values together, as far as they do not need a value that needs the
// branch -- another branch of the same condition then follows that value. The two sides of a
// condition follow each other, which a backend makes an `if` and its `else`.
//
// Of the values that can be computed, the first in a depth-first order from the roots in their
// order comes first: every root right after what it needs, so that a value is computed close to
// where it is used. With the writes of the lowering, which run as soon as their values are computed
// (Store::position), the values of roots that need little of each other -- the peaks of RotD50
// along 180 angles, each the maximum of its old value and the velocity along its angle -- live one
// after the other, not all at once: the live values, which a GPU keeps in its registers and spills
// beyond them, are those that roots share (the velocities) and few more.
//
// A backend that runs a branch only for the points that take it computes the same values as one
// that runs it for all, and its loads are as safe: every offset of a load is clamped into its array
// for every coordinate, NaN and the infinities included (Expand.cpp), so a load the hardware
// executes for a lane that does not take the branch -- predicated, speculated or masked off --
// stays in its array. No instruction of a stage communicates between points (there is no shuffle,
// vote, barrier or shared memory; a check reports with a plain store of one value), so the points
// of a GPU warp may diverge at a branch and converge after it in any order.

#include "Expr/Ir.h"
#include "Expr/Program.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <tuple>
#include <vector>

namespace seissol::expr {

/// The points of the region `parent` where `condition` holds (`side` true) or does not; region 0,
/// every point, has none.
struct Region {
  std::int32_t parent{-1};
  NodeId condition{NoNode};
  bool side{true};
  std::int32_t depth{0};
};

/// The region of every value of a program, or of a stage of it (see REGIONS above).
class Regions {
  public:
  /// The regions of `values`, which are in ascending NodeId order and hold the operands of each of
  /// them but of those `leaves` marks (by NodeId, if given), which have none here; `roots` use
  /// their values in region 0. Without `truths`, a select whose condition is a truth value
  /// (ValueType::Bool) uses its operands where it lies, as any other operation.
  Regions(const Arena& arena,
          const std::vector<NodeId>& values,
          const std::vector<NodeId>& roots,
          const std::vector<char>* leaves = nullptr,
          bool truths = true);

  /// The region of `value`; -1 for a value no root needs.
  [[nodiscard]] std::int32_t of(NodeId value) const;
  [[nodiscard]] const Region& operator[](std::int32_t region) const {
    return regions_[static_cast<std::size_t>(region)];
  }
  [[nodiscard]] std::int32_t size() const { return static_cast<std::int32_t>(regions_.size()); }
  /// The innermost region that holds both.
  [[nodiscard]] std::int32_t common(std::int32_t x, std::int32_t y) const;
  /// The region in which `user`, of the region `region`, uses its operand number `operand` (in the
  /// order of Arena::children()): `region`, or a branch of it for the operands of a select.
  [[nodiscard]] std::int32_t
      useOf(const Arena& arena, NodeId user, std::int32_t region, std::size_t operand) const;

  private:
  std::int32_t branch(std::int32_t parent, NodeId condition, bool side);
  void meet(NodeId value, std::int32_t region);

  [[nodiscard]] bool branches(const Arena& arena, const Node& node) const;

  bool truths_{true};
  std::vector<Region> regions_;
  std::map<std::tuple<std::int32_t, NodeId, bool>, std::int32_t> branches_;
  std::vector<std::int32_t> regionOf_;
};

/// Where the outputs, the states and the checks of a program use its values, exactly (see WHERE A
/// VALUE IS USED above).
class Reach {
  public:
  /// The uses of the values of `program`, for truth values built in `target` (another arena) of
  /// the type `type`, ValueType::Bool, or Real in a program not yet typed: `condition` gives the
  /// condition of a select of `program` in `target`.
  Reach(const Program& program,
        Arena& target,
        std::function<NodeId(NodeId)> condition,
        ValueType type);

  /// A value of `target` that is `value` at the points where the program uses `node`, and
  /// `otherwise` at the others: selects over the conditions of the regions around `node`, so that
  /// `value` lies in its region, and over those that tell the points of the region apart where
  /// `node` is not used.
  [[nodiscard]] NodeId masked(NodeId node, NodeId value, NodeId otherwise);

  /// Whether the program uses `node` at every point, where masked() is `value` itself; without
  /// building a value, so before the conditions are.
  [[nodiscard]] bool everywhere(NodeId node);

  private:
  struct Use {
    NodeId user{NoNode}; // NoNode: a root
    std::int32_t region{0};
  };

  /// The points of the region of `node` where it is used, as a truth value of `target`; NoNode for
  /// all of them.
  NodeId within(NodeId node);
  NodeId withinFromUses(NodeId node);
  /// `x` where `condition` is `side`, else `y`
  NodeId choose(NodeId condition, bool side, NodeId x, NodeId y);

  const Arena& arena_;
  Arena& target_;
  std::function<NodeId(NodeId)> condition_;
  ValueType type_;
  std::vector<NodeId> roots_;
  std::vector<NodeId> values_;
  Regions regions_;
  std::vector<std::vector<Use>> uses_;
  std::vector<char> done_;
  std::vector<NodeId> within_;
  // whether within_ is NoNode, where known: -1 not yet
  std::vector<signed char> plain_;
  NodeId true_{NoNode};
  NodeId false_{NoNode};
};

/// One step of the code of a stage: a value; or the start of a branch, of the points where the
/// condition `node` is `side`; or its end.
struct Step {
  enum class Kind : std::uint8_t { Value, If, EndIf };
  Kind kind{Kind::Value};
  NodeId node{NoNode};
  bool side{true};
};

/// `values`, as Regions takes them with `roots` and `leaves`, in an order to compute them in, with
/// a branch around each region whose values weigh `threshold` or more by `weight` (see THE ORDER
/// above); in the depth-first order from `roots` where no region does.
[[nodiscard]] std::vector<Step> branchedOrder(const Arena& arena,
                                              const std::vector<NodeId>& values,
                                              const std::vector<NodeId>& roots,
                                              const std::vector<char>& leaves,
                                              const std::function<std::int64_t(NodeId)>& weight,
                                              std::int64_t threshold);

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_BRANCHES_H_
