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
// NAMING: this is what the previous discussion called "field"/"FieldStore". The
// word is taken: `Kind::Field` in the IR is a named input channel, inherited
// from the derived-output port where "field" means vx/syy/…. Two different
// things called Field in adjacent headers is a bug waiting to happen, so the
// external-data concept is "grid" throughout. Say the word if you'd rather keep
// FieldStore and rename the IR node instead — it is a one-line change either
// way, but it has to be decided before packages 1–3 start.
//
// SCOPE: fully replicated, host memory, one copy per rank (or per node via the
// existing shared-memory tier in AsagiLite). No device mirroring, no sharding,
// no texture objects. Time-dependent grids hold a sliding window of slices.

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace seissol::datafield {

enum class GridKind : std::uint8_t { AsagiLite, Scec, Distributed /* not implemented */ };

enum class Interpolation : std::uint8_t { Nearest, Linear, Cubic };

// What happens outside the grid volume. Clamp is the default because it is the
// only option that lets the sub-box/halo work later: with Fail, a point one
// stencil width outside the local box is an error rather than a well-defined
// edge value.
enum class BoundaryMode : std::uint8_t { Clamp, Fail };

// Stencil half-width per interpolation kind: 1, 2, 4 points per axis.
constexpr std::size_t stencilWidth(Interpolation interp) {
  switch (interp) {
  case Interpolation::Nearest:
    return 1;
  case Interpolation::Linear:
    return 2;
  case Interpolation::Cubic:
    return 4;
  }
  return 1;
}

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
  // This is what the Lookup node lowers to. The stencil gather inside is a
  // scalar loop (random access), the tensor-product reduction over it is not —
  // implementations should keep those two phases separate so the second one
  // vectorises.
  virtual void sampleBatch(const double* x, std::size_t count, double* out) const = 0;

  // Volume, for the future sub-box analysis and for diagnostics.
  virtual void bounds(double* min, double* max) const = 0;

  // Time window currently resident. Absent for static grids.
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
  // be called before any resident window expires:
  //   min over time-dependent grids of (stencilWidth - 1) * dt_file.
  // Returns nullopt when no grid is time-dependent — in which case the caller
  // should not register the update module at all.
  [[nodiscard]] std::optional<double> suggestedSyncInterval() const;

  private:
  std::vector<GridDesc> descs_;
  std::vector<std::unique_ptr<Grid>> grids_;
};

} // namespace seissol::datafield

#endif // SEISSOL_SRC_READER_DATAFIELD_GRID_H_
