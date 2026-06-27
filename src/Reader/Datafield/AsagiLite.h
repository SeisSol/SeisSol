// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_DATAFIELD_ASAGILITE_H_
#define SEISSOL_SRC_READER_DATAFIELD_ASAGILITE_H_

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <mpi.h>
#include <string>
#include <type_traits>

namespace asagi_lite {

enum class Error {
  Success = 0,
  IoError,
  BadRank,
  BadElementType,
  AllocationFailed,
};

/// Rounding mode used when converting a continuous coordinate to a grid index.
enum class Rounding : std::uint8_t {
  Down,    // floor
  Nearest, // round half away from zero
  Up,      // ceil
};

/// Drop-in replacement for an ASAGI grid that loads a NetCDF-4/HDF5 file
/// once into memory and serves local lookups afterwards.
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
///                          threshold `kSharedMemThreshold` is used.
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
class Grid {
  public:
  /// Up to 3 spatial dimensions plus one component dimension.
  static constexpr std::size_t KMaxFileDimensions = 4;
  /// Hard fallback threshold used by the auto-tier logic when node memory
  /// cannot be queried.
  static constexpr std::size_t KSharedMemThreshold = 16ULL << 30; // 16 GiB

  enum class ElementType : std::uint8_t { Float, Double };
  enum class Owner : std::uint8_t { None, Malloc, Shmem };

  Grid() noexcept = default;
  /// Collective on the node communicator in shared-memory mode.
  ~Grid() noexcept;

  Grid(const Grid&) = delete;
  Grid& operator=(const Grid&) = delete;
  Grid(Grid&& other) noexcept;
  Grid& operator=(Grid&& other) noexcept;

  /// Defaults to `MPI_COMM_SELF`. For cross-rank sharing set e.g.
  /// `MPI_COMM_WORLD` before calling `open()`.
  void setComm(MPI_Comm comm) noexcept { comm_ = comm; }

  Error open(const std::string& filename, const std::string& varname = "z");
  void close() noexcept;
  [[nodiscard]] bool isOpen() const noexcept { return data_ != nullptr; }

  [[nodiscard]] std::size_t dimensions() const noexcept { return ndims_; }
  [[nodiscard]] std::size_t sizeAt(std::size_t i) const noexcept;
  [[nodiscard]] double minAt(std::size_t i) const noexcept;
  [[nodiscard]] double maxAt(std::size_t i) const noexcept;
  [[nodiscard]] ElementType elementType() const noexcept { return elemType_; }
  [[nodiscard]] Owner owner() const noexcept { return owner_; }

  // ---------- Coordinates -> indices ----------

  template <Rounding R = Rounding::Nearest>
  [[nodiscard]] std::size_t coordToIndex(std::size_t i, double coord) const noexcept;

  template <std::size_t N>
  std::array<std::size_t, N> coordsToIndices(const std::array<double, N>& coords) const noexcept;

  template <std::size_t N>
  std::array<std::size_t, N>
      coordsToIndices(const std::array<double, N>& coords,
                      const std::array<Rounding, N>& rounding) const noexcept;

  /// Pointer variant for runtime-determined dimensionality.
  /// Pass `rounding == nullptr` to use Rounding::Nearest for all dimensions.
  void coordsToIndices(const double* coords,
                       const Rounding* rounding,
                       std::size_t n,
                       std::size_t* outIndices) const noexcept;

  // ---------- Indices -> values ----------

  template <typename T, std::size_t N, std::size_t M = 1>
  std::array<T, M> getValue(const std::array<std::size_t, N>& indices) const noexcept;

  template <typename T>
  void getValue(const std::size_t* indices,
                std::size_t nIndices,
                T* out,
                std::size_t mComponents) const noexcept;

  private:
  Error allocateLocal(std::size_t bytes);
  Error allocateShared(std::size_t bytes);
  void freeBuffer() noexcept;
  void resetMetadata() noexcept;

  [[nodiscard]] double coordPos(std::size_t i, double coord) const noexcept {
    return (coord - offset_[i]) * invScaling_[i];
  }

  std::array<std::size_t, KMaxFileDimensions> size_{{1, 1, 1, 1}};
  std::array<std::size_t, KMaxFileDimensions> stride_{{0, 0, 0, 0}};
  std::array<double, KMaxFileDimensions> offset_{{0.0, 0.0, 0.0, 0.0}};
  std::array<double, KMaxFileDimensions> invScaling_{{1.0, 1.0, 1.0, 1.0}};

  std::size_t ndims_ = 0;
  void* data_ = nullptr;
  ElementType elemType_ = ElementType::Float;
  Owner owner_ = Owner::None;

  MPI_Comm comm_ = MPI_COMM_SELF;
  MPI_Win win_ = MPI_WIN_NULL;
  MPI_Comm nodeComm_ = MPI_COMM_NULL;   // ranks sharing one node
  MPI_Comm leaderComm_ = MPI_COMM_NULL; // one rank per node (broadcast)
};

// ============================ Hot path ============================

template <Rounding R>
inline std::size_t Grid::coordToIndex(std::size_t i, double coord) const noexcept {
  assert(i < ndims_ && "dimension index out of range");
  const double pos = coordPos(i, coord);
  if constexpr (R == Rounding::Down) {
    return static_cast<std::size_t>(std::floor(pos));
  } else if constexpr (R == Rounding::Up) {
    return static_cast<std::size_t>(std::ceil(pos));
  } else {
    return static_cast<std::size_t>(std::lround(pos));
  }
}

template <std::size_t N>
inline std::array<std::size_t, N>
    Grid::coordsToIndices(const std::array<double, N>& coords) const noexcept {
  static_assert(N >= 1 && N <= KMaxFileDimensions, "N must be between 1 and kMaxFileDimensions");
  std::array<std::size_t, N> out{};
  for (std::size_t i = 0; i < N; ++i) {
    out[i] = coordToIndex<Rounding::Nearest>(i, coords[i]);
  }
  return out;
}

template <std::size_t N>
inline std::array<std::size_t, N>
    Grid::coordsToIndices(const std::array<double, N>& coords,
                          const std::array<Rounding, N>& rounding) const noexcept {
  static_assert(N >= 1 && N <= KMaxFileDimensions, "N must be between 1 and kMaxFileDimensions");
  std::array<std::size_t, N> out{};
  for (std::size_t i = 0; i < N; ++i) {
    const double pos = coordPos(i, coords[i]);
    switch (rounding[i]) {
    case Rounding::Down:
      out[i] = static_cast<std::size_t>(std::floor(pos));
      break;
    case Rounding::Up:
      out[i] = static_cast<std::size_t>(std::ceil(pos));
      break;
    case Rounding::Nearest:
    default:
      out[i] = static_cast<std::size_t>(std::lround(pos));
      break;
    }
  }
  return out;
}

template <typename T, std::size_t N, std::size_t M>
inline std::array<T, M> Grid::getValue(const std::array<std::size_t, N>& indices) const noexcept {
  static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "T must be float or double");
  static_assert(N >= 1 && N <= KMaxFileDimensions, "N must be between 1 and kMaxFileDimensions");
  static_assert(M >= 1, "M must be at least 1");

  assert(data_ != nullptr && "grid is not open");
  assert(((M == 1 && ndims_ == N) || (ndims_ == N + 1 && size_[N] == M)) &&
         "grid shape does not match getValue<T,N,M> parameters");

  std::size_t idx = 0;
  for (std::size_t i = 0; i < N; ++i) {
    assert(indices[i] < size_[i] && "index out of range");
    idx += indices[i] * stride_[i];
  }

  std::array<T, M> out;
  if (elemType_ == ElementType::Float) {
    const float* p = static_cast<const float*>(data_) + idx;
    for (std::size_t j = 0; j < M; ++j) {
      out[j] = static_cast<T>(p[j]);
    }
  } else {
    const double* p = static_cast<const double*>(data_) + idx;
    for (std::size_t j = 0; j < M; ++j) {
      out[j] = static_cast<T>(p[j]);
    }
  }
  return out;
}

extern template void
    Grid::getValue<float>(const std::size_t*, std::size_t, float*, std::size_t) const noexcept;
extern template void
    Grid::getValue<double>(const std::size_t*, std::size_t, double*, std::size_t) const noexcept;

} // namespace asagi_lite
#endif // SEISSOL_SRC_READER_DATAFIELD_ASAGILITE_H_
