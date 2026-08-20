// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_DATAFIELD_GRID_H_
#define SEISSOL_SRC_READER_DATAFIELD_GRID_H_

// External data grids (ASAGI-style structured grids) for the expression layer.
//
// NAMING: this is what earlier drafts called "field"/"FieldStore". The word is
// taken: `Kind::Field` in the IR is a named input channel, inherited from the
// derived-output port where "field" means vx/syy/…. Two different things called
// Field in adjacent headers is a bug waiting to happen, so the external-data
// concept is "grid" throughout.
//
// SCOPE: fully replicated, host memory, one copy per rank (or per node via the
// shared-memory tier in AsagiLite). No device mirroring, no sharding, no texture
// objects. Time-dependent grids hold a sliding window of slices.
//
// The interpolation model is the one from easi's InterpolationKernel.h
// (davschneller/precommit): every scheme is separable, so the d-dimensional
// interpolant is the tensor product of d one-dimensional interpolants, and a
// scheme is fully described by its 1-D stencil geometry plus its 1-D weight
// function. Porting that file rather than reinventing it keeps the two code
// paths numerically identical, which matters because the same .nc file can be
// read through easi and through here in the same run.

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace seissol::datafield {

enum class GridKind : std::uint8_t { AsagiLite, Scec, Distributed /* not implemented */ };

// Nearest: O(h), reproduces constants.
// Linear:  O(h^2), C^0, reproduces degree-1 polynomials.
// Cubic:   O(h^3), C^1, reproduces degree-2 polynomials. Keys' cubic
//          convolution with a = -1/2 (Catmull-Rom).
//
// WARNING: Cubic is not bound-preserving. Across a sharp contrast it overshoots
// by roughly 7% of the jump height. For a velocity model with a Moho-type
// discontinuity that can produce non-physical rho or mu — a negative or
// near-zero value will not be caught here, it will surface much later as a
// nonsensical wave speed or a diverging time step. Use Linear or Nearest
// wherever the interpolated values must stay within the sampled range, and
// treat Cubic as opt-in per grid rather than as a global default.
enum class Interpolation : std::uint8_t { Nearest, Linear, Cubic };

// Geometry of a one-dimensional stencil. With
//   base = floor((x - min) / delta)
// the stencil covers grid indices [base + offset, base + offset + width - 1].
// If roundToNearest is set, base is obtained by rounding and the stencil
// degenerates to that single point.
struct StencilKernel {
  unsigned width{2};
  int offset{0};
  bool roundToNearest{false};
};

constexpr StencilKernel kernelOf(Interpolation interp) {
  switch (interp) {
  case Interpolation::Nearest:
    return {1, 0, true};
  case Interpolation::Linear:
    return {2, 0, false};
  case Interpolation::Cubic:
    return {4, -1, false};
  }
  return {2, 0, false};
}

inline constexpr unsigned MaxStencilWidth = 4;

// What happens when a query point lies OUTSIDE the grid volume.
//
// This is not the same question as "the stencil does not fit near the edge".
// That case needs no mode: inside the volume the scheme always degrades to the
// linear fallback with zero weights on the slots that fall outside, so nothing
// is ever extrapolated. BoundaryMode only decides what an out-of-volume query
// does, which today is undefined behaviour in AsagiLite — coordToIndex casts a
// negative std::floor result to std::size_t and reads far out of bounds.
//
// Clamp is the default because it is the only option compatible with a future
// sub-box/halo decomposition, where a point one stencil width outside the local
// box must still yield a well-defined edge value rather than an error.
enum class BoundaryMode : std::uint8_t { Clamp, Fail };

struct GridDesc {
  std::string path;
  std::string variable;
  GridKind kind{GridKind::AsagiLite};
  Interpolation interpolation{Interpolation::Linear};
  BoundaryMode boundary{BoundaryMode::Clamp};

  // Which file axis is time, if any. Absent = static grid, loaded once.
  std::optional<int> timeAxis;

  [[nodiscard]] bool timeDependent() const { return timeAxis.has_value(); }
  [[nodiscard]] bool operator==(const GridDesc& other) const;
};

// One replicated grid. Coordinates are always f64 at this interface regardless
// of the program's compute type; the conversion happens in the Lookup lowering,
// because index arithmetic in f32 loses exactness well before the value does.
class Grid {
  public:
  Grid() = default;
  virtual ~Grid() = default;
  Grid(const Grid&) = delete;
  Grid& operator=(const Grid&) = delete;
  Grid(Grid&&) = delete;
  Grid& operator=(Grid&&) = delete;

  [[nodiscard]] virtual std::size_t dimensions() const = 0;
  [[nodiscard]] virtual std::size_t components() const = 0;

  // Single point. `x` has dimensions() entries, `out` has components().
  virtual void sample(const double* x, double* out) const = 0;

  // Batch, SoA in and SoA out: x[axis * count + lane], out[comp * count + lane].
  // This is what a Lookup node lowers to.
  //
  // Implementations must keep the two phases separate: gathering the width^d
  // stencil values is random access and will not vectorise, whereas the
  // tensor-product reduction over them will. Fusing them costs the reduction
  // its vectorisation, which is where most of the work is for Cubic (4^3 = 64
  // values per point in 3-D).
  virtual void sampleBatch(const double* x, std::size_t count, double* out) const = 0;

  // Volume, for the out-of-bounds decision and for diagnostics.
  virtual void bounds(double* min, double* max) const = 0;

  // Time window currently resident, as [first, last] sample times. Absent for
  // static grids.
  [[nodiscard]] virtual std::optional<std::pair<double, double>> timeWindow() const {
    return std::nullopt;
  }
};

// Owns every grid a set of Programs references, deduplicated by GridDesc.
class GridStore {
  public:
  GridStore() = default;
  ~GridStore() = default;
  GridStore(const GridStore&) = delete;
  GridStore& operator=(const GridStore&) = delete;
  GridStore(GridStore&&) = delete;
  GridStore& operator=(GridStore&&) = delete;

  // Returns a stable index; opening the file is deferred to load().
  [[nodiscard]] std::size_t intern(const GridDesc& desc);
  [[nodiscard]] const Grid& get(std::size_t id) const;
  [[nodiscard]] std::size_t size() const { return grids_.size(); }

  // How many time slices a time-dependent grid keeps resident. Must exceed the
  // widest stencil in use, see suggestedSyncInterval().
  void setResidentSlices(std::size_t slices);
  [[nodiscard]] std::size_t residentSlices() const { return residentSlices_; }

  // Opens every interned grid. Collective over the rank's communicator.
  void load();

  // Advance every time-dependent grid so that its resident window covers
  // `time`. Idempotent and cheap when nothing needs reloading, which is what
  // lets the caller just invoke it at every synchronisation point.
  //
  // There is deliberately no user-facing sampling interval. The cadence is a
  // property of the file, not of the configuration: the time axis has its own
  // spacing, and the window stays valid until `time` leaves it. Asking the user
  // to also configure an interval would only create a second, independent way
  // to get it wrong.
  void update(double time);

  // Largest synchronisation interval at which update() is still guaranteed to
  // run before any resident window expires.
  //
  // With W resident slices and a stencil {w, o}, a query whose base index is b
  // (relative to the window start) is servable exactly when
  //     b + o >= 0   and   b + o + w - 1 <= W - 1,
  // i.e. for b in [-o, W - w - o] — that is W - w + 1 cells, spanning
  // (W - w + 1) * dtFile of simulated time. One cell of headroom is kept so the
  // reload happens strictly before expiry rather than exactly at it, giving
  //     min over time-dependent grids of (W - w) * dtFile.
  //
  // Returns nullopt when no grid is time-dependent — in which case the caller
  // should not register the update module at all. Throws std::invalid_argument
  // if W <= w for any grid, since then no interval would be safe.
  [[nodiscard]] std::optional<double> suggestedSyncInterval() const;

  private:
  std::vector<GridDesc> descs_;
  std::vector<std::unique_ptr<Grid>> grids_;
  std::size_t residentSlices_{MaxStencilWidth + 1};
};

} // namespace seissol::datafield

#endif // SEISSOL_SRC_READER_DATAFIELD_GRID_H_
