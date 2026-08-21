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
// shared-memory tier in AsagiLiteGrid). No device mirroring, no sharding, no texture
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

namespace seissol::reader::datafield {

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

// Matches easi's Grid<Derived>::MaxDimensions. Do not raise one without the
// other; the two paths are compared value-for-value (see the tolerance note).
inline constexpr unsigned MaxGridDimensions = 6;

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
//
// The clamp must be written so that a non-finite coordinate lands on the
// boundary rather than overflowing the cast, i.e. `if (!(raw > 0.0)) raw = 0.0;`
// and not `if (raw < 0.0)`. NaN fails the first test and passes the second.
enum class BoundaryMode : std::uint8_t { Clamp, Fail };

// ============================================================================
// AXIS ORDER  --  decided, do not re-derive per backend
// ============================================================================
//
// Axis 0 is the FASTEST-VARYING file dimension, i.e. easi's / ASAGI's
// convention, and the reverse of the HDF5 shape order. For a variable declared
// `z(y, x)` the grid axes are (x, y).
//
// AsagiLite stores the file order directly (`stride_[ndims_-1] == 1`, axis 0
// slowest) and is therefore the side that reverses. That reversal happens once,
// in the AsagiLiteGrid adapter, not at every call site. The component axis (see
// `components()`) is stripped BEFORE the reversal, so it never appears as a
// spatial axis.
//
// This is the single most dangerous thing in this header: reading a velocity
// model transposed produces a plausible-looking field and a wrong simulation.
// The cross-check against easi is what proves it, and it is not optional.

// ============================================================================
// GEOMETRY DERIVATION  --  decided, do not re-derive per backend
// ============================================================================
//
// `delta` is not stored in a NetCDF file; it is derived from the coordinate
// variable, and the two obvious derivations do not agree:
//
//   first spacing      delta = v[1] - v[0]                  (easi / ASAGI)
//   endpoint average   delta = (v[n-1] - v[0]) / (n - 1)    (AsagiLite today)
//
// Coordinate variables are routinely stored as float32, and on such an axis the
// two differ by far more than a rounding. Measured on the nominal node
// positions of a float32 axis:
//
//   lat 45.25 deg, 15 arcsec, 4801 pts : max offset 1.17 cells, 2618/4801 nodes
//                                        land in a different base cell
//   lon 5.50 deg, 30 arcsec, 2401 pts  : max offset 0.037 cells, 1120/2401
//   depth 0, 0.1 km, 5001 pts          : max offset 7.4e-5 cells, 1380/5001
//   x -300 km, 250 m, 4001 pts         : exact in binary, 0/4001
//
// A whole-cell offset means reading the wrong material. DECIDED: the easi form
// wins. `delta` is the first spacing and the reciprocal is formed once as
// `1.0 / delta`, not as `(n - 1) / span` -- the latter differs from the former
// by one ulp on about a quarter of all inputs.
//
// `num` stays the HDF5 extent rather than easi's `lround((max-min)/delta) + 1`,
// because the extent is known exactly here and the reconstruction is not. Where
// the two disagree the coordinate axis is not uniform to within half a spacing,
// the uniform-grid model does not hold for that file, and the backend should
// logWarning rather than silently pick one.

// ============================================================================
// CROSS-PATH TOLERANCE  --  why this is not a bitwise criterion
// ============================================================================
//
// easi is a separate library with its own build flags. With FMA contraction
// enabled on one side and not the other, the tensor-product reduction -- a sum
// of width^d weighted products -- associates differently and the results differ.
//
// Measured, gcc 13 -O2 -march=x86-64-v3, -ffp-contract=off vs =fast, 8192
// values per (dimension, scheme):
//
//   anchored on the RESULT : up to 84089 ulp, max relative 1.4e-11
//   anchored on the DATA   : |delta| <= 4.3 * eps * max|v| , in every case
//
// The first number is not a defect, it is the wrong metric: the interpolant has
// zero crossings, and near one of them any relative or ulp measure of an
// absolute perturbation diverges. The perturbation itself scales with the
// STENCIL VALUES, not with the result.
//
// The cross-path acceptance test therefore compares against
//
//   tol = GridCrossPathTolerance * eps * max|v over the stencil|
//
// The analytic worst case is width^d * sum|w| = 64 * 1.25^3 ~ 125 for Cubic in
// 3-D; 16 leaves an order of magnitude over what is observed while staying well
// inside the analytic bound, so a real divergence still trips it.
inline constexpr double GridCrossPathTolerance = 16.0;

// Default total memory for the resident time windows of all time-dependent
// grids on one rank. Deliberately modest: a run that wants long windows should
// say so, and a run with no time-dependent grids pays nothing either way.
inline constexpr std::size_t DefaultWindowMemoryBudget = 1ULL << 30; // 1 GiB

// Uniform Cartesian geometry of one grid, in grid axis order (see above).
//
// Everything the scheduling side needs is derivable from this: the volume is
// [min, min + (num-1)*delta], and the file's time step is delta[timeAxis]. That
// is why there is no separate bounds()/timeSpacing() virtual -- one export, one
// place for a backend to get it wrong.
struct GridGeometry {
  unsigned dimensions{0};
  double min[MaxGridDimensions]{};
  double delta[MaxGridDimensions]{};
  unsigned num[MaxGridDimensions]{};
};

struct GridDesc {
  std::string path;
  std::string variable;
  GridKind kind{GridKind::AsagiLite};
  Interpolation interpolation{Interpolation::Linear};
  BoundaryMode boundary{BoundaryMode::Clamp};

  // Which GRID axis is time, if any -- grid numbering, not file numbering.
  // Absent = static grid, loaded once.
  std::optional<int> timeAxis;

  [[nodiscard]] bool timeDependent() const { return timeAxis.has_value(); }
  [[nodiscard]] bool operator==(const GridDesc& other) const;
};

// NOT PROVIDED, deliberately: a component-name to component-index mapping.
// easi carries one (ASAGI::permutation) because its component names come from
// the model file. Here `Kind::Lookup` already carries a resolved integer
// component, and resolving a name against a grid is the caller's business --
// the script or the parameter file says which component it means. Adding a
// second, independent name-resolution path here is exactly the class of
// "the name silently resolved to the wrong column" defect that Program.h
// documents as the reason the signature is inferred rather than declared.

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

  // Uniform geometry in grid axis order. Backends fill dimensions() entries.
  virtual void geometry(GridGeometry& out) const = 0;

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
  //
  // The reduction must NOT run in place. easi's scalar driver does, and it is
  // safe there only because the accumulator is a register: the k == 0 source
  // slot and the destination slot are the same, and the write happens after all
  // reads. Vectorised over lanes there is no such register, so an in-place pass
  // zeroes the k == 0 term before reading it. Ping-pong between two buffers.
  // The failure is invisible in 1-D (a single pass never aliases) and wrong in
  // every value in 2-D and 3-D.
  virtual void sampleBatch(const double* x, std::size_t count, double* out) const = 0;

  // Volume, for the out-of-bounds decision and for diagnostics. Derived from
  // geometry(), not virtual -- a backend that could disagree with its own
  // geometry is a backend that will.
  void bounds(double* min, double* max) const;

  // Bytes one slice of the time axis costs in memory — the whole array divided
  // by the number of time slices. Feeds GridStore's window sizing. Zero means
  // the backend cannot say, in which case the store falls back to the smallest
  // window the scheme allows.
  [[nodiscard]] virtual std::size_t bytesPerTimeSlice() const { return 0; }

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

  // How many time slices each time-dependent grid keeps resident.
  //
  // By default this is DERIVED from a memory budget rather than configured: the
  // budget is split evenly over the time-dependent grids and each one takes as
  // many slices as its share buys, floored at the narrowest window its scheme
  // admits (width + 1) and capped at its own number of slices. That turns a
  // magic constant into a physical quantity — the user has an opinion about how
  // much memory a velocity model may occupy, not about how many slices that is.
  //
  // Sizing per grid rather than globally matters because the cost of a slice
  // varies by orders of magnitude between grids: a 2-D surface field and a 3-D
  // volume field in the same run should not be forced to the same window.
  //
  // setResidentSlices() overrides the derivation for every grid, for the cases
  // where someone does have an opinion.
  void setWindowMemoryBudget(std::size_t bytes);
  [[nodiscard]] std::size_t windowMemoryBudget() const { return windowMemoryBudget_; }
  void setResidentSlices(std::optional<std::size_t> slices);
  [[nodiscard]] std::optional<std::size_t> residentSlicesOverride() const {
    return residentSlicesOverride_;
  }
  /// Slices grid `id` keeps resident: the override if set, else the budget
  /// derivation. Requires load() to have run, since it reads the geometry.
  [[nodiscard]] std::size_t residentSlicesFor(std::size_t id) const;

  // Fallback used when a grid's time axis carries no usable spacing -- no
  // coordinate variable, a single slice, or a non-finite delta. Suboptimal by
  // construction: it is a global number standing in for a per-file property,
  // so a grid that falls back gets a logWarning naming the file. Comes from the
  // parameter file, not from getenv.
  void setDefaultTimeSpacing(double dt);
  [[nodiscard]] double defaultTimeSpacing() const { return defaultTimeSpacing_; }

  // The file's time step for grid `id`: geometry().delta[*desc.timeAxis], or
  // defaultTimeSpacing() where that is not usable. Absent for static grids.
  [[nodiscard]] std::optional<double> timeSpacing(std::size_t id) const;

  // Opens every interned grid. Collective over the rank's communicator.
  void load();

  /// Installs a backend directly, bypassing load(). Exists so that the store's
  /// arithmetic — deduplication, window sizing, the sync interval — can be
  /// tested without a file on disk, since none of it depends on where the bytes
  /// came from. Production code calls load().
  void injectForTesting(std::size_t id, std::unique_ptr<Grid> grid);

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
  // Both inputs come from what is already here: dtFile is timeSpacing(id), and
  // w is kernelOf(descs_[id].interpolation).width -- the interpolation scheme is
  // per grid, not per axis, so the time axis uses the same stencil as space.
  // (If nearest-in-time with cubic-in-space is ever wanted, that is a second
  // Interpolation field on GridDesc, not a change here.)
  //
  // W here is residentSlicesFor(id), so the interval follows the memory budget:
  // more memory buys a longer window buys fewer synchronisation points. The
  // caller is expected to take the minimum of this and whatever its own
  // scheduling requires — this number is an upper bound on a safe interval, not
  // a request.
  //
  // Returns nullopt when no grid is time-dependent — in which case the caller
  // should not register the update module at all. Throws std::invalid_argument
  // if W <= w for any grid, since then no interval would be safe.
  [[nodiscard]] std::optional<double> suggestedSyncInterval() const;

  private:
  std::vector<GridDesc> descs_;
  std::vector<std::unique_ptr<Grid>> grids_;
  std::optional<std::size_t> residentSlicesOverride_;
  std::size_t windowMemoryBudget_{DefaultWindowMemoryBudget};
  double defaultTimeSpacing_{0.0};
};

} // namespace seissol::reader::datafield

#endif // SEISSOL_SRC_READER_DATAFIELD_GRID_H_
