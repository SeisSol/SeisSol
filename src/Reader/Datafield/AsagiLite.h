// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_DATAFIELD_ASAGILITE_H_
#define SEISSOL_SRC_READER_DATAFIELD_ASAGILITE_H_

#include <array>
#include <cstddef>
#include <cstdint>
// This header owns raw HDF5 handles, so it needs hdf5.h. That is acceptable
// because nothing outside the factory includes it: consumers go through Grid.
#include "Grid.h"
#include "Interpolation.h"

#include <hdf5.h>
#include <mpi.h>
#include <optional>
#include <string>
#include <utility>

namespace seissol::reader::datafield {

/// Replacement for an ASAGI grid that loads a NetCDF-4/HDF5 file once into
/// memory and serves local lookups afterwards.
///
/// MPI is a hard requirement; consumers must link against MPI (and typically
/// against parallel HDF5).
///
/// Storage tier is selected automatically at `open()` time. Two environment
/// variables override the decision:
///   - `ASAGI_LITE_USE_SHMEM`
///       unset (default) -> auto: use MPI-3 shared memory iff replicating
///                          the dataset across every rank on the node would
///                          exceed roughly half of physical node memory.
///                          If memory cannot be queried, the hard fallback
///                          threshold `SharedMemThreshold` is used.
///       `=1`/`on`/...   -> force shared memory.
///       `=0`/`off`/...  -> force per-rank malloc.
///   - `ASAGI_LITE_BCAST`
///       unset (default) -> every reader rank reads the file itself.
///       `=1`/`on`/...   -> global rank 0 reads; the result is broadcast.
///                          In shared-memory mode the broadcast targets a
///                          leader-only sub-communicator (one rank per node);
///                          otherwise it targets the full communicator.
///
/// Supported file layouts (HDF5/NetCDF dimension order, slowest first):
///   - Scalar grid: shape `[d_0, ..., d_{N-1}]`, M = 1
///   - Vector grid: shape `[d_0, ..., d_{N-1}, M]`
///
/// AXIS ORDER: the file order above is the reverse of the grid axis order this
/// class exposes through Grid (see Grid.h). The reversal happens exactly once,
/// in buildView(), by permuting stride/num/min/deltaInv into `view_`. Nothing on
/// the sampling path reverses anything, so there is no per-call opportunity to
/// get it half right.
///
/// TIME WINDOWING: a time-dependent grid keeps only `residentSlices` slices of
/// its time axis in memory and slides that window forward as the simulation
/// advances. The file stays open for the lifetime of the grid.
///
/// The time axis must be the SLOWEST file dimension — grid axis
/// dimensions() - 1. That is the NetCDF convention (`z(t, y, x)`) and it is what
/// makes a window a single contiguous hyperslab, so sliding it is one read and
/// one memmove rather than a strided gather. A file with time anywhere else is
/// rejected with a message saying so rather than silently read at a fraction of
/// the speed.
///
/// Sliding reuses the overlap: the slices the old and new window share are
/// memmoved down and only the genuinely new ones are read. The alternative, a
/// ring buffer, would avoid the memmove but put modular arithmetic on the time
/// axis into the innermost sampling loop, which is a bad trade for a window of
/// a handful of slices.
class AsagiLiteGrid final : public Grid {
  public:
  /// Up to 3 spatial dimensions plus one component dimension.
  static constexpr std::size_t MaxFileDimensions = 4;
  /// Hard fallback threshold used by the auto-tier logic when node memory
  /// cannot be queried.
  static constexpr std::size_t SharedMemThreshold = 16ULL << 30; // 16 GiB

  enum class Owner : std::uint8_t { None, Malloc, Shmem };

  /// Failure modes of open(). Internal: the Grid interface has no error return,
  /// so makeGrid() turns anything but Success into a logError naming the file.
  enum class Error : std::uint8_t {
    Success = 0,
    IoError,
    BadRank,
    BadElementType,
    AllocationFailed,
  };

  AsagiLiteGrid() noexcept = default;
  /// Collective on the node communicator in shared-memory mode.
  ~AsagiLiteGrid() noexcept override;

  AsagiLiteGrid(const AsagiLiteGrid&) = delete;
  AsagiLiteGrid& operator=(const AsagiLiteGrid&) = delete;
  AsagiLiteGrid(AsagiLiteGrid&&) = delete;
  AsagiLiteGrid& operator=(AsagiLiteGrid&&) = delete;

  /// Defaults to `MPI_COMM_SELF`. For cross-rank sharing set e.g.
  /// `MPI_COMM_WORLD` before calling `open()`.
  void setComm(MPI_Comm comm) noexcept { comm_ = comm; }

  /// `interp` and `boundary` come from the GridDesc and are fixed for the
  /// lifetime of the grid; sampleBatch() has no way to take them per call and
  /// should not, since the scheme is a property of the model, not of the query.
  /// Opens the file and reads its metadata. A static grid is fully loaded here.
  /// A time-dependent one is not: call resizeWindow() once the window size is
  /// known, which is what GridStore::load() does after applying its budget.
  Error open(const std::string& filename,
             const std::string& varname = "z",
             Interpolation interp = Interpolation::Linear,
             BoundaryMode boundary = BoundaryMode::Clamp,
             std::optional<int> timeAxis = std::nullopt);

  /// Allocates the resident window and loads its initial slices. Collective.
  /// No-op for a static grid, which open() already loaded.
  Error resizeWindow(std::size_t residentSlices);

  /// Slides the window so that `time` is servable at the full accuracy of the
  /// configured scheme. Collective. Returns without touching the file when the
  /// current window already serves it, which is the common case.
  Error advanceTo(double time);

  [[nodiscard]] TimeWindow window() const noexcept { return window_; }
  void close() noexcept;
  [[nodiscard]] bool isOpen() const noexcept { return data_ != nullptr; }

  [[nodiscard]] Owner owner() const noexcept { return owner_; }
  [[nodiscard]] ElementType elementType() const noexcept { return elemType_; }

  // ---------- Grid interface ----------

  [[nodiscard]] std::size_t dimensions() const override { return view_.dimensions; }
  [[nodiscard]] std::size_t components() const override { return view_.components; }
  void geometry(GridGeometry& out) const override;
  void sample(const double* x, double* out) const override;
  void sampleBatch(const double* const* x, std::size_t count, double* const* out) const override;
  void sampleBatch(const double* const* x, std::size_t count, float* const* out) const override;

  private:
  /// Shared by both sampleBatch overloads: the BoundaryMode::Fail diagnostic is
  /// about the query, not about the result type.
  void reportClamped(std::size_t clamped, std::size_t count) const;

  public:
  [[nodiscard]] std::size_t bytesPerTimeSlice() const override { return sliceBytes_; }
  [[nodiscard]] std::optional<std::pair<double, double>> timeWindow() const override;

  private:
  Error allocateLocal(std::size_t bytes);
  Error allocateShared(std::size_t bytes);
  void freeBuffer() noexcept;
  void resetMetadata() noexcept;
  /// Permutes the file-order metadata into `view_` and strips the component
  /// axis. Called once at the end of open().
  void buildView() noexcept;
  /// Reads `count` file slices starting at `first` into slice slot `destSlot`.
  /// Collective; honours the same reader/broadcast tiering as the initial load.
  [[nodiscard]] Error readSlices(std::size_t first, std::size_t count, std::size_t destSlot);
  /// Element index of the first entry of a resident slice.
  [[nodiscard]] std::size_t sliceOffset(std::size_t slot) const noexcept {
    return slot * sliceElems_;
  }

  // ---------- file-order metadata, slowest dimension first ----------
  std::array<std::size_t, MaxFileDimensions> size_{{1, 1, 1, 1}};
  std::array<std::size_t, MaxFileDimensions> stride_{{0, 0, 0, 0}};
  std::array<double, MaxFileDimensions> offset_{{0.0, 0.0, 0.0, 0.0}};
  /// Reciprocal of the first spacing, formed as 1.0 / (v[1] - v[0]). See the
  /// geometry-derivation note in Grid.h for why it is not (n-1)/(v[n-1]-v[0]).
  std::array<double, MaxFileDimensions> deltaInv_{{1.0, 1.0, 1.0, 1.0}};
  std::array<double, MaxFileDimensions> delta_{{1.0, 1.0, 1.0, 1.0}};
  /// Set when the last file dimension carries no HDF5 dimension scale, which is
  /// how a vector grid is told apart from an N+1-dimensional scalar grid.
  bool hasComponentAxis_ = false;

  std::size_t ndims_ = 0;
  /// Kept open for the lifetime of the grid: a sliding window has to read again.
  hid_t fileId_ = -1;
  hid_t dsetId_ = -1;
  /// Resident portion of the time axis. count == numTimeSlices_ for a static
  /// grid or one that fits entirely, in which case advanceTo() never reads.
  TimeWindow window_;
  std::size_t numTimeSlices_ = 1;
  std::size_t sliceElems_ = 0;
  std::size_t sliceBytes_ = 0;
  /// Origin of the time axis in the FILE, as opposed to view_.min[timeAxis],
  /// which tracks the window and therefore moves.
  double timeOrigin_ = 0.0;
  void* data_ = nullptr;
  ElementType elemType_ = ElementType::Float;
  Owner owner_ = Owner::None;

  // ---------- grid-order view handed to the sampler ----------
  ArrayView view_;
  /// Batch working memory, reused across calls. Mutable because sampleBatch()
  /// is const on the Grid interface and must stay so -- the scratch is not part
  /// of the observable state.
  mutable SampleScratch scratch_;
  Interpolation interp_ = Interpolation::Linear;
  BoundaryMode boundary_ = BoundaryMode::Clamp;
  std::optional<int> timeAxis_;
  std::string path_;

  MPI_Comm comm_ = MPI_COMM_SELF;
  MPI_Win win_ = MPI_WIN_NULL;
  MPI_Comm nodeComm_ = MPI_COMM_NULL;   // ranks sharing one node
  MPI_Comm leaderComm_ = MPI_COMM_NULL; // one rank per node (broadcast)
};

} // namespace seissol::reader::datafield
#endif // SEISSOL_SRC_READER_DATAFIELD_ASAGILITE_H_
