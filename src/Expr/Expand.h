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
// backend and the same as those of datafield::sampleBatch(). The whole numbers -- the indices, the
// offsets -- are exact in binary64.
//
// The geometry is read when the kernel runs, not baked into it: the window of a time-dependent grid
// moves it, and a kernel samples whichever file of the same scheme it is bound to.
//
// A grid whose boundary mode fails beyond its volume gets a check per lookup (GridCheck): a value
// that is nonzero where a coordinate was clamped onto the volume, which every backend reports.

#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Interpolation.h"

#include <cstdint>
#include <vector>

namespace seissol::expr {

/// A value that is nonzero where a lookup of `grid`, a grid that fails beyond its volume, samples
/// beyond it.
struct GridCheck {
  GridId grid{NoGrid};
  NodeId root{NoNode};
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

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_EXPAND_H_
