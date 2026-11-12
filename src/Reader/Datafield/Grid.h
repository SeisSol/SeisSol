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
// shared-memory tier in AsagiLiteGrid). No sharding, no texture objects.
// Time-dependent grids hold a sliding window of slices. A copy elsewhere -- the
// one a device kernel samples -- is a GridStore::Mirror of the array a grid
// exposes (Grid::arrayView()), which the store keeps current as the window moves.
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

// Inline: a grid whose values and knots the program gives itself (GridData), e.g. the layers of a
// layered model.
enum class GridKind : std::uint8_t { AsagiLite, Scec, Distributed /* not implemented */, Inline };

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
//
// Previous and Next: the steps of a layered model (easi's `lower` and `upper`), the value of the
// last node below the coordinate, or of the first node at or above it -- a coordinate on a node
// belongs to the cell below it, as in easi. Beyond the axis the nearest end. Only the written-out
// lookup (expr::expandLookups()) samples them, as it samples an axis of knots.
enum class Interpolation : std::uint8_t { Nearest, Linear, Cubic, Previous, Next };

// The nodes of an axis: uniform (its first node, its spacing and their number), or knots, a list
// of strictly increasing, finite coordinates -- a layered model, a file of non-uniform coordinates.
// A lookup locates a coordinate among knots by a binary search.
enum class AxisKind : std::uint8_t { Uniform, Knots };

// The most knots an axis of knots of a file can have: what the binary search of a written-out
// lookup covers, whose steps are fixed when the program is lowered, before the file is read.
inline constexpr std::size_t MaxKnotsOfFile = std::size_t{1} << 24;

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
  case Interpolation::Previous:
  case Interpolation::Next:
    return {1, 0, true};
  }
  return {2, 0, false};
}

inline constexpr unsigned MaxStencilWidth = 4;

// Slices a window of a time-dependent grid keeps to spare on either side of the stencils of the
// times it is placed for, those up to the next update (placeTimeWindow(),
// GridStore::suggestedSyncInterval()). Those times never come near the ends of the window, so a
// written-out lookup checks that the window holds its time with the predicates its stencil
// computes anyway, which are not exact at the ends (expr::CheckKind::Window): no binary64
// operation for the check, and no tolerance.
inline constexpr unsigned TimeWindowSlack = 1;

// The narrowest window of a time-dependent grid with the stencil `kernel`: the stencil, the slices
// to spare on either side, and one cell to move on, so the window is reloaded at every update but
// never too late.
constexpr std::size_t narrowestWindow(StencilKernel kernel) {
  return std::size_t{kernel.width} + 1 + 2 * std::size_t{TimeWindowSlack};
}

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
//   first spacing      delta = v[1] - v[0]
//   endpoint average   delta = (v[n-1] - v[0]) / (n - 1)    (ASAGI, and so easi)
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
// A whole-cell offset means reading the wrong material. DECIDED: the form easi
// samples with wins. easi takes the spacing of ASAGI (NetCdfReader::getScaling(),
// (last - first) / (size - 1)) and of its own regular grid ((max - min) / (num -
// 1)): the endpoint average, which is also the one whose error does not grow
// along the axis. The reciprocal is formed once as `1.0 / delta`, not as
// `(n - 1) / span` -- the latter differs from the former by one ulp on about a
// quarter of all inputs.
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
// Surface-syntax names for the two closed vocabularies above. Both frontends
// (sderiv grid declarations, Lua field_specs) spell these as strings; a mapping
// in each would let "linear" and Interpolation drift apart between the traced
// and the interpreted path. One mapping, here, next to the enums it maps onto.
//
// std::nullopt rather than a default on an unknown name: the sderiv frontend
// reports it with a source position, the Lua tracer as a TraceFailure, and
// neither diagnostic can be produced from in here.
[[nodiscard]] std::optional<GridKind> parseGridKind(const std::string& text);
[[nodiscard]] std::optional<Interpolation> parseInterpolation(const std::string& text);
/// "clamp" or "fail".
[[nodiscard]] std::optional<BoundaryMode> parseBoundaryMode(const std::string& text);
/// "uniform" or "knots".
[[nodiscard]] std::optional<AxisKind> parseAxisKind(const std::string& text);
/// The spellings of the vocabularies: what the parsers above read.
[[nodiscard]] const char* spelling(GridKind kind);
[[nodiscard]] const char* spelling(Interpolation interpolation);
[[nodiscard]] const char* spelling(BoundaryMode boundary);
[[nodiscard]] const char* spelling(AxisKind axis);

// The dataset inside the file when a frontend does not name one. Matches
// AsagiLiteGrid::open's own default; stated once so the two frontends cannot
// disagree about what an omitted variable means.
inline constexpr const char* DefaultGridVariable = "z";

// The array of a grid, which a compiled kernel samples (Interpolation.h).
struct ArrayView;

struct GridGeometry {
  unsigned dimensions{0};
  double min[MaxGridDimensions]{};
  double delta[MaxGridDimensions]{};
  unsigned num[MaxGridDimensions]{};
};

/// What the file of a grid says about one of its axes, where the backend reads one: the name of
/// its dimension, and the attributes units, positive and standard_name -- empty where absent.
struct GridAxisInfo {
  std::string dimension;
  std::string units;
  std::string positive;
  std::string standardName;
};

/// The values and the knots of a grid that a program gives itself (GridKind::Inline).
struct GridData {
  /// per axis of the grid its knots, strictly increasing and finite
  std::vector<std::vector<double>> knots;
  /// the values at the nodes: the components of a node together, the nodes in the order of the
  /// axes, the first fastest
  std::vector<double> values;
  [[nodiscard]] bool operator==(const GridData& other) const;
};

struct GridDesc {
  std::string path;
  std::string variable;
  GridKind kind{GridKind::AsagiLite};
  /// how every axis interpolates and what a sample beyond it does, unless `interpolations` and
  /// `boundaries` say so per axis
  Interpolation interpolation{Interpolation::Linear};
  BoundaryMode boundary{BoundaryMode::Clamp};
  /// per axis, as many as the grid has axes, or none: the same for all
  std::vector<Interpolation> interpolations;
  std::vector<BoundaryMode> boundaries;
  /// per axis whether it is uniform or knots (of a file: from its coordinate variable); none for
  /// all uniform. Every axis of an inline grid is knots.
  std::vector<AxisKind> axes;
  /// the values and the knots of an inline grid
  std::shared_ptr<const GridData> data;

  [[nodiscard]] Interpolation interpolationOf(std::size_t axis) const {
    return axis < interpolations.size() ? interpolations[axis] : interpolation;
  }
  [[nodiscard]] BoundaryMode boundaryOf(std::size_t axis) const {
    return axis < boundaries.size() ? boundaries[axis] : boundary;
  }
  [[nodiscard]] AxisKind axisOf(std::size_t axis) const {
    if (kind == GridKind::Inline) {
      return AxisKind::Knots;
    }
    return axis < axes.size() ? axes[axis] : AxisKind::Uniform;
  }
  /// Whether a sample beyond the volume fails on some axis.
  [[nodiscard]] bool failsAnywhere() const;

  // Which GRID axis is time, if any -- grid numbering, not file numbering.
  // Absent = static grid, loaded once.
  std::optional<int> timeAxis;

  // The names a script gave the grid and its components, which the core format of sderiv prints
  // (Program::gridNames()). Descriptive only: neither the comparison of two grids nor the key of a
  // kernel has them, so that they travel with every copy of the grid without changing what it is.
  std::string name;
  std::vector<std::string> componentNames;

  [[nodiscard]] bool timeDependent() const { return timeAxis.has_value(); }
  [[nodiscard]] bool operator==(const GridDesc& other) const;
  /// What a message calls the grid: its file, or the name of an inline grid, which has none.
  [[nodiscard]] std::string label() const { return path.empty() ? name : path; }

  // Everything about this grid that reaches the EMITTED KERNEL, and nothing
  // else. Program::canonicalForm() folds this in so that two programs sampling
  // different files through the same scheme share a compiled kernel, while two
  // sampling the same file through Linear and Cubic do not.
  //
  // path and variable are deliberately absent: they select which bytes the
  // GridStore loads, not which code the backend emits. Program::internGrid()
  // therefore dedupes on canonicalKey() PLUS path PLUS variable, and the cache
  // key on canonicalKey() alone -- two different questions, two different keys.
  //
  // Stable across runs and builds: fixed-width tokens, no enum values written
  // as integers (a reordered enum must not silently alias two schemes), no
  // locale-dependent formatting.
  //
  // The values and the knots of an inline grid are in it: they are the program, as its constants
  // are, and two tables of other values are two programs.
  [[nodiscard]] std::string canonicalKey() const;
};

/// A grid as a frontend declares it, before it is checked: the strings and the numbers of its
/// attributes (the sderiv specification, `grid`). The interpolation and the boundary are one for
/// all axes or one per axis; `axes` gives the kind of each axis of a file; `knots` and `values`
/// are the data of an inline grid.
struct GridSpec {
  std::string kind;
  std::string path;
  std::string variable;
  std::vector<std::string> components;
  std::vector<std::string> interpolation;
  std::vector<std::string> boundary;
  std::vector<std::string> axes;
  std::vector<std::vector<double>> knots;
  std::vector<double> values;
  std::optional<int> timeAxis;
};

/// The GridDesc of `spec`, its components named as it names them. Throws std::invalid_argument,
/// with what is wrong, for a spec that does not describe a grid.
[[nodiscard]] GridDesc gridDescOf(const GridSpec& spec);

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

  // Batch. One pointer per axis in, one per component out, each addressing a
  // contiguous run of `count`: x[axis][lane], out[comp][lane]. `out` has
  // components() entries; a null entry skips that component entirely.
  //
  // A lookup of a program does not come here: the lowering writes it out
  // (Expr/Expand.h) into the arithmetic below, on the array of arrayView(), so
  // that every backend computes it with its own instructions -- and two lookups
  // of one grid at the same coordinates share the cells and the weights. This is
  // the reference that arithmetic is checked against bit for bit, and what the
  // other consumers of a grid call.
  //
  // Pointer arrays rather than flat blocks, because a caller typically holds each
  // coordinate in its own run and wants each result in its own, so flat blocks
  // would cost a copy of dimension * count values in and components * count out,
  // purely to satisfy the signature.
  //
  // The null-skip is the other half of the same decision: a caller that wants one
  // component would otherwise pay the full tensor-product reduction for values
  // nobody reads, 4^3 weighted adds per point per component for Cubic in 3-D.
  // Because componentStride is 1, the random-access address arithmetic is shared
  // across components, which is why this is one call and not one per component.
  //
  // ALIASING: out[comp] may alias x[axis], because the slot allocator can reuse
  // a coordinate slot that dies at the Lookup as its destination. Implementations
  // must therefore read every coordinate before writing any output. Distinct
  // out[] entries must not alias each other.
  //
  // The float form differs only in the final store; the reduction stays f64 in
  // both, because the cross-path tolerance above is derived against easi's f64
  // kernel and an f32 reduction lands far outside it.
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
  virtual void sampleBatch(const double* const* x, std::size_t count, double* const* out) const = 0;
  virtual void sampleBatch(const double* const* x, std::size_t count, float* const* out) const = 0;

  // Volume, for the out-of-bounds decision and for diagnostics. Derived from
  // geometry(), not virtual -- a backend that could disagree with its own
  // geometry is a backend that will.
  void bounds(double* min, double* max) const;

  // Bytes one slice of the time axis costs in memory — the whole array divided
  // by the number of time slices. Feeds GridStore's window sizing. Zero means
  // the backend cannot say, in which case the store falls back to the smallest
  // window the scheme allows.
  [[nodiscard]] virtual std::size_t bytesPerTimeSlice() const { return 0; }

  // What the file says about the axes of the grid, in grid axis order; empty where the backend
  // reads no file.
  [[nodiscard]] virtual std::vector<GridAxisInfo> axisInfo() const { return {}; }

  // The slices of the file along the time axis, of which the window holds some
  // (geometry() describes the window). Zero where the backend cannot say, in
  // which case the store takes the geometry for the whole file.
  [[nodiscard]] virtual std::size_t timeSlices() const { return 0; }

  // Time window currently resident, as [first, last] sample times. Absent for
  // static grids.
  [[nodiscard]] virtual std::optional<std::pair<double, double>> timeWindow() const {
    return std::nullopt;
  }

  // The array the samples come from, where the backend keeps one in memory: what a compiled kernel
  // samples itself, by the arithmetic of Interpolation.h, and what a device copies. Null where the
  // backend keeps none. The address stays the same for the lifetime of the grid; what it holds
  // changes where revision() does.
  [[nodiscard]] virtual const ArrayView* arrayView() const { return nullptr; }

  // Changes whenever the array of arrayView() or its geometry does: a window that slides or is
  // resized. A copy of the array is current as long as this has not changed since it was taken.
  [[nodiscard]] virtual std::uint64_t revision() const { return 0; }
};

// Owns every grid a set of Programs references, deduplicated by GridDesc.
class GridStore;

// The process-wide store.
//
// Injectable so tests stay hermetic: pass an own GridStore to CompiledReader
// rather than reaching for this.
[[nodiscard]] GridStore& sharedGridStore();

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
  /// The SHA-256 of the file of grid `id`, as load() read it on the first rank (fileHashOf());
  /// empty before, for a grid without a file, and where hashing is off.
  [[nodiscard]] const std::string& fileHash(std::size_t id) const { return hashes_.at(id); }
  [[nodiscard]] std::size_t size() const { return grids_.size(); }

  // How many time slices each time-dependent grid keeps resident.
  //
  // By default this is DERIVED from a memory budget rather than configured: the
  // budget is split evenly over the time-dependent grids and each one takes as
  // many slices as its share buys, floored at the narrowest window its scheme
  // admits (narrowestWindow()) and capped at its own number of slices. That turns a
  // magic constant into a physical quantity — the user has an opinion about how
  // much memory a velocity model may occupy, not about how many slices that is.
  //
  // Sizing per grid rather than globally matters because the cost of a slice
  // varies by orders of magnitude between grids: a 2-D surface field and a 3-D
  // volume field in the same run should not be forced to the same window.
  //
  // setResidentSlices() overrides the derivation for every grid, for the cases
  // where someone does have an opinion. A window narrower than the stencil of a
  // grid is refused by suggestedSyncInterval(), per grid.
  void setWindowMemoryBudget(std::size_t bytes);
  [[nodiscard]] std::size_t windowMemoryBudget() const { return windowMemoryBudget_; }
  void setResidentSlices(std::optional<std::size_t> slices);
  [[nodiscard]] std::optional<std::size_t> residentSlicesOverride() const {
    return residentSlicesOverride_;
  }
  /// Slices grid `id` keeps resident: the override if set, else the budget
  /// derivation; never more than its file has. Requires load() to have run, since
  /// it reads the geometry.
  [[nodiscard]] std::size_t residentSlicesFor(std::size_t id) const;
  /// The slices of the file of grid `id` along its time axis (1 for a grid without one).
  [[nodiscard]] std::size_t timeSlicesOf(std::size_t id) const;

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

  // Opens every interned grid that is not open yet, and sizes the window of every time-dependent
  // grid whose resident slices changed (residentSlicesFor()). Collective over the rank's
  // communicator wherever it opens or sizes anything, and only there: a call that finds nothing to
  // do -- the call of every kernel made after the grids it samples were loaded -- touches neither
  // the files nor the communicator. That is what lets a kernel be made on some ranks only (those
  // with a face of a boundary, say) once every rank has interned and loaded the grids of all
  // scripts in the same order, which the initialization does as it loads the scripts.
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
  //
  // The mirrors follow the grids that changed (Mirror::refresh()).
  void update(double time);

  // A copy of a grid elsewhere -- in the memory of a device, say -- that has to follow the grid as
  // its window moves. The store owns it, and refreshes it after every load() and update(), which
  // are the only places a grid changes.
  class Mirror {
public:
    Mirror() = default;
    virtual ~Mirror() = default;
    Mirror(const Mirror&) = delete;
    Mirror& operator=(const Mirror&) = delete;
    Mirror(Mirror&&) = delete;
    Mirror& operator=(Mirror&&) = delete;

    /// The grid may have changed: compare Grid::revision() with that of the copy.
    virtual void refresh(const Grid& grid) = 0;
    /// The end of the run (finish()): what the copy has to report.
    virtual void finish(const Grid& /*grid*/) {}
  };

  /// The mirror of grid `id` of the kind `key` (an own number per kind of copy), or null.
  [[nodiscard]] Mirror* mirror(std::size_t id, int key) const;
  /// Attaches `mirror` to grid `id` as its copy of the kind `key`, in place of one there.
  Mirror& attachMirror(std::size_t id, int key, std::unique_ptr<Mirror> mirror);

  /// The end of the run: every mirror reports what it has to (Mirror::finish()), and the grids
  /// close their files, which a store that lives until the process ends could not do later, past
  /// HDF5. They keep their arrays and sample on -- the analysis of the initial condition after the
  /// end of the simulation evaluates its script, say -- but their windows do not move any more.
  void finish();

  // Largest synchronisation interval at which update() is still guaranteed to
  // run before any resident window expires.
  //
  // With W resident slices and a stencil {w, o}, a query whose base index is b
  // (relative to the window start) is servable exactly when
  //     b + o >= 0   and   b + o + w - 1 <= W - 1.
  // A window keeps s = TimeWindowSlack slices to spare on either side, so the
  // queries until the next update have b in [s - o, W - w - o - s] -- that is
  // W - w - 2s + 1 cells, spanning (W - w - 2s + 1) * dtFile of simulated time.
  // One cell of headroom is kept so the reload happens strictly before expiry
  // rather than exactly at it, giving
  //     min over time-dependent grids of (W - w - 2s) * dtFile.
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
  // A grid whose window holds its whole file never moves and asks for no
  // interval. Returns nullopt when no grid is time-dependent, or none moves —
  // in which case the caller should not register the update module at all.
  // Throws std::invalid_argument if W <= w + 2s for any grid, since then no
  // interval would be safe.
  [[nodiscard]] std::optional<double> suggestedSyncInterval() const;

  private:
  void refreshMirrors();

  std::vector<GridDesc> descs_;
  std::vector<std::unique_ptr<Grid>> grids_;
  std::vector<std::string> hashes_;
  /// The resident slices each time-dependent grid was sized to by load(); zero before.
  std::vector<std::size_t> windowSlices_;
  /// Per grid, its mirrors and their kinds. After grids_, so that they go first.
  std::vector<std::vector<std::pair<int, std::unique_ptr<Mirror>>>> mirrors_;
  std::optional<std::size_t> residentSlicesOverride_;
  std::size_t windowMemoryBudget_{DefaultWindowMemoryBudget};
  double defaultTimeSpacing_{0.0};
};

} // namespace seissol::reader::datafield

#endif // SEISSOL_SRC_READER_DATAFIELD_GRID_H_
