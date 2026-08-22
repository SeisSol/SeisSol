// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_DATAFIELD_INTERPOLATION_H_
#define SEISSOL_SRC_READER_DATAFIELD_INTERPOLATION_H_

// Separable interpolation on a replicated uniform Cartesian array.
//
// Ported from easi's InterpolationKernel.h / component/Grid.h
// (davschneller/precommit). The weight functions and the locate branches are
// transcribed unchanged: the same .nc file can be read through easi and through
// here in the same run, and the two results are compared against
// GridCrossPathTolerance. Any edit here is an edit to that comparison.
//
// Deliberately free of HDF5 and MPI so the numerics can be unit-tested without
// a file. Everything specific to the backend — where the bytes came from, which
// file dimension maps to which grid axis — is baked into the ArrayView the
// caller hands over.

#include "Grid.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace seissol::reader::datafield {

enum class ElementType : std::uint8_t { Float, Double };

/// One-dimensional weights of the scheme at local coordinate s in [0, 1].
/// Writes kernelOf(type).width entries.
void interpolationWeights(Interpolation type, double s, double* weights);

/// Weights that make a width-w stencil act like linear interpolation between
/// slots `slotOfBase` and `slotOfBase + 1`. Used where the full stencil does not
/// fit inside the array, i.e. within (width/2 - 1) cells of either edge. This is
/// what keeps the scheme from extrapolating there; the order drops locally, and
/// that is the intended behaviour, not a fallback to be optimised away.
void linearFallbackWeights(unsigned width, unsigned slotOfBase, double s, double* weights);

/// A replicated, strided array with uniform Cartesian geometry, in GRID axis
/// order (axis 0 = fastest-varying file dimension, see Grid.h).
///
/// The axis permutation lives entirely in `stride` and `num`: a backend whose
/// file order differs from grid order permutes these once at open() time. The
/// sampling loops below never reverse anything, which is the only way to keep
/// the transposition bug from having somewhere to hide.
struct ArrayView {
  const void* data{nullptr};
  ElementType type{ElementType::Float};
  unsigned dimensions{0};
  unsigned components{1};

  double min[MaxGridDimensions]{};
  /// Reciprocal of the axis spacing, formed once as 1.0 / delta — see the
  /// geometry-derivation note in Grid.h. Zero marks a degenerate axis
  /// (num == 1, or a coordinate variable that does not increase), which then
  /// always samples index 0.
  double deltaInv[MaxGridDimensions]{};
  unsigned num[MaxGridDimensions]{};
  /// Element stride along each grid axis.
  std::size_t stride[MaxGridDimensions]{};
  /// Element stride between consecutive components of one grid point.
  std::size_t componentStride{1};
};

/// Per-batch working memory. Owned by the caller so that a sampling loop does
/// not allocate; grows monotonically and is reused.
class SampleScratch {
  public:
  void reserve(unsigned dimensions, unsigned width, unsigned components, std::size_t count);

  std::vector<double> weights;     // [(axis * width + slot) * count + lane]
  std::vector<std::int64_t> start; // [axis * count + lane]
  std::vector<double> gather;      // [(flat * components + c) * count + lane]
  std::vector<double> work;        // ping-pong partner of `gather`
};

/// Batch sample. One pointer per axis and one per component, each addressing a
/// contiguous run of `count` values: x[axis][lane], out[component][lane].
///
/// CHANGED (reported): both used to be flat blocks, x[axis * count + lane] and
/// out[component * count + lane]. Pointer arrays because every caller that
/// matters holds these runs SEPARATELY -- the interpreter keeps each coordinate
/// in its own transient slot and wants each result in its own -- so a flat block
/// forced a copy of dimension * count values in and components * count out, per
/// tile, purely to satisfy the signature.
///
/// `out[c] == nullptr` skips component c entirely: it is neither gathered nor
/// reduced. This is not a micro-optimisation. A Kind::Lookup node names ONE
/// component, so without skipping, every lookup would pay the full
/// tensor-product reduction for components nobody reads -- 4^3 weighted adds
/// per point per component for Cubic in 3-D. Gathering all components at a
/// stencil point stays worthwhile in the other direction, which is why the
/// split is here and not one call per component: componentStride is 1, so the
/// random-access address arithmetic is shared, and it is the expensive half.
///
/// ALIASING: out[c] may alias x[d]. The interpreter's slot allocator is free to
/// reuse a coordinate slot that dies at this instruction as the destination.
/// This is safe because every read of x happens in phase 1, before the first
/// write to out -- a constraint on the implementation, not an accident.
/// out[c] pointers must not alias EACH OTHER; that is a caller error.
///
/// Coordinates outside the array volume are clamped onto the boundary; the
/// return value counts how many lanes were clamped on at least one axis, which
/// is what a BoundaryMode::Fail backend turns into a logError. Non-finite
/// coordinates are clamped too rather than being allowed to reach the index
/// cast (see the note on the clamp idiom in Grid.h).
///
/// The float form differs ONLY in the final store. The reduction itself stays
/// f64 in both: the cross-path tolerance in Grid.h is derived against easi's
/// f64 kernel, and reducing in f32 would move the result far outside it.
std::size_t sampleBatch(const ArrayView& view,
                        Interpolation interp,
                        const double* const* x,
                        std::size_t count,
                        double* const* out,
                        SampleScratch& scratch);
std::size_t sampleBatch(const ArrayView& view,
                        Interpolation interp,
                        const double* const* x,
                        std::size_t count,
                        float* const* out,
                        SampleScratch& scratch);

/// Single point. Convenience wrapper; the batch form is the one that vectorises.
std::size_t sample(const ArrayView& view, Interpolation interp, const double* x, double* out);

// ===========================================================================
// Time window placement
// ===========================================================================
//
// A time-dependent grid keeps `count` consecutive slices of its time axis
// resident, starting at file slice `start`. Kept here rather than in the
// backend because it is pure arithmetic on the stencil geometry, and because
// getting it wrong degrades accuracy silently: a mispositioned window still
// produces a number, just one interpolated by the edge fallback instead of the
// requested scheme.

struct TimeWindow {
  std::size_t start{0};
  std::size_t count{0};

  [[nodiscard]] bool operator==(const TimeWindow& o) const {
    return start == o.start && count == o.count;
  }
};

/// Whether `base` — a cell index on the file's time axis — is inside `window` at
/// the FULL accuracy of `kernel`, i.e. without falling back to linear.
///
/// Near either end of the file the fallback is unavoidable and correct, exactly
/// as it is in space; those positions report true so the window does not chase
/// a placement that does not exist.
bool windowServes(const TimeWindow& window,
                  std::int64_t base,
                  StencilKernel kernel,
                  std::size_t numSlices);

/// Where the window has to sit to serve `base`.
///
/// The window is placed so that `base` lands on the FIRST slot the scheme can
/// use, putting the whole remaining width ahead of the query. Simulated time
/// only ever moves forward, so a centred window would throw away half its
/// headroom and reload twice as often. The usable span is `count - width` cells,
/// which is where suggestedSyncInterval()'s (W - w) * dt comes from.
TimeWindow placeTimeWindow(std::int64_t base,
                           StencilKernel kernel,
                           std::size_t resident,
                           std::size_t numSlices);

} // namespace seissol::reader::datafield

#endif // SEISSOL_SRC_READER_DATAFIELD_INTERPOLATION_H_
