// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_EXPAND_H_
#define SEISSOL_SRC_EXPR_EXPAND_H_

// The lookups of a program, written out (part of the lowering).
//
// A frontend samples a grid with one node, Kind::Lookup. The lowering writes it out as what it
// computes: per axis, the cell and the weights of the stencil of the scheme from the coordinate and
// the geometry of the grid (Kind::GridParam); the elements of the stencil (Kind::GridLoad); and
// their tensor product. Every backend then computes a sample from the same instructions as any
// other value -- the interpreter, the compiled CPU kernels and the device kernels alike, with no
// sampler of their own -- and the lowering sees into it: two lookups of one grid at the same
// coordinates share the cells and the weights (e.g. rho, vp and vs of a material), and so do the
// lookups at the same position at several times; a part that does not change between calls is
// hoisted like any other.
//
// THE ARITHMETIC IS THAT OF Interpolation.cpp, operation by operation and in the same order (the
// branches of locateAxis as selects, which give the same values), in binary64 whatever the compute
// type: the coordinates are converted to binary64 first, and the sample to the type of the lookup
// last. Its multiplications, additions and subtractions are strict (Node::strict), so that a device
// kernel does not contract them into fused multiply-adds: the samples are the same bits on every
// backend and the same as those of datafield::sampleBatch().
//
// The indices, the node counts and the offsets are whole numbers (ValueType::I64), the conditions
// truth values (ValueType::Bool), as in Interpolation.cpp: exact, and on a device with little
// binary64 throughput far cheaper than binary64. The offsets leave out the component, which the
// load adds (Kind::GridLoad), so that the lookups of several components at one point share them.
//
// The geometry is read when the kernel runs, not baked into it: the window of a time-dependent grid
// moves it, and a kernel samples whichever file of the same scheme it is bound to.
//
// A grid whose boundary mode fails beyond its volume gets a check per lookup (GridCheck): true
// where a coordinate was clamped onto the volume, which every backend reports. A time-dependent
// grid gets one on its time axis: true where the window does not hold the stencil of the time,
// where the file goes on beyond it -- a time the window was not placed for.

#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Interpolation.h"

#include <cstdint>
#include <string>
#include <vector>

namespace seissol::expr {

/// What a check of a written-out lookup finds where it is true.
enum class CheckKind : std::uint8_t {
  /// a sample beyond the volume of a grid whose boundary mode is to fail
  Outside,
  /// a sample at a time beyond the resident window of a time-dependent grid, where the file goes
  /// on: a script that samples the grid at other times than those from one synchronization point
  /// to the next, for which the window is placed, or an update of the window that did not run
  Window,
  /// a point that does not meet a check of the program (Program::checks()): its condition is zero
  User
};

/// A truth value (Bool) that is true where a lookup of `grid` finds what `kind` says, or, of kind
/// User, where a point does not meet the check `index` of the program.
struct GridCheck {
  GridId grid{NoGrid};
  CheckKind kind{CheckKind::Outside};
  NodeId root{NoNode};
  std::int32_t index{-1};
};

/// How a program samples a grid: with how many coordinates (0 for a grid it does not sample), and
/// up to which component. A backend checks the grid it is bound to against it.
struct GridUse {
  std::int32_t dimensions{0};
  std::int32_t components{0};
};

struct ExpandedProgram {
  /// the program without lookups
  Program program;
  std::vector<GridCheck> checks;
  /// by GridId
  std::vector<GridUse> uses;
};

/// `program`, typed (typed()), with its lookups written out. Throws std::invalid_argument for a
/// grid sampled with different numbers of coordinates.
[[nodiscard]] ExpandedProgram expandLookups(const Program& program);

/// How `program` samples its grids (by GridId), as expandLookups() finds it, without writing the
/// lookups out. Throws std::invalid_argument for a grid sampled with different numbers of
/// coordinates.
[[nodiscard]] std::vector<GridUse> gridUses(const Program& program);

/// The value a Kind::GridParam reads of `view`.
[[nodiscard]] double
    gridField(const reader::datafield::ArrayView& view, GridField field, std::int32_t axis);

/// Raises what a check of `kind` found for the lookups of the grid at `path` (logError()), which
/// every backend reports alike.
void reportCheck(CheckKind kind, const std::string& path);

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_EXPAND_H_
