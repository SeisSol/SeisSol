// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "AsagiLite.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <mpi.h>
#include <optional>
#include <string>
#include <string_view>
#include <sys/types.h>
#include <unistd.h> // sysconf

namespace asagi_lite {
namespace {

// ---------- Environment ----------

/// Tri-state parse: unset -> nullopt, "0"/"off"/... -> false, otherwise true.
std::optional<bool> envBool(const char* name) noexcept {
  const char* v = std::getenv(name);
  if ((v == nullptr) || ((*v) == 0)) {
    return std::nullopt;
  }
  const std::string_view s{v};
  if (s == "0" || s == "false" || s == "FALSE" || s == "no" || s == "off") {
    return false;
  }
  return true;
}

bool envOn(const char* name) noexcept { return envBool(name).value_or(false); }

// ---------- Node memory probing ----------

/// Returns ~half of physical node memory as a budget, or 0 if unavailable.
/// On non-Linux systems where `_SC_PHYS_PAGES` is missing this falls back to
/// returning 0, in which case the caller uses `kSharedMemThreshold`.
std::size_t nodeMemoryBudget() noexcept {
#if defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
  const long pages = sysconf(_SC_PHYS_PAGES);
  const long pageSize = sysconf(_SC_PAGESIZE);
  if (pages > 0 && pageSize > 0) {
    return (static_cast<std::size_t>(pages) * static_cast<std::size_t>(pageSize)) / 2;
  }
#endif
  return 0;
}

/// Auto-tier decision.  Splits `comm` into a node-local sub-communicator and
/// hands it back to the caller via `outNodeComm`; the caller reuses it if
/// shmem is chosen, or frees it otherwise.
bool decideShmem(std::size_t bytes, MPI_Comm comm, MPI_Comm& outNodeComm) noexcept {
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &outNodeComm);

  if (const auto override = envBool("ASAGI_LITE_USE_SHMEM"); override.has_value()) {
    return *override;
  }

  int ranksOnNode = 1;
  MPI_Comm_size(outNodeComm, &ranksOnNode);
  const std::size_t replicated = bytes * static_cast<std::size_t>(ranksOnNode);

  if (const std::size_t budget = nodeMemoryBudget(); budget > 0) {
    return replicated > budget;
  }
  // No memory information: fall back to the fixed per-rank threshold.
  return bytes >= Grid::KSharedMemThreshold;
}

// ---------- MPI helpers ----------

/// MPI_Bcast count is `int`; chunk to handle buffers larger than ~2 GiB.
void bcastLarge(void* buf, std::size_t bytes, int root, MPI_Comm comm) noexcept {
  constexpr std::size_t KChunk = std::size_t{1} << 30; // 1 GiB
  auto* p = static_cast<char*>(buf);
  while (bytes > 0) {
    const int n = static_cast<int>(std::min(bytes, KChunk));
    MPI_Bcast(p, n, MPI_BYTE, root, comm);
    p += n;
    bytes -= static_cast<std::size_t>(n);
  }
}

// ---------- HDF5 RAII ----------

template <herr_t (*Closer)(hid_t)>
class HidGuard {
  public:
  explicit HidGuard(hid_t id) noexcept : id_(id) {}
  ~HidGuard() noexcept {
    if (id_ >= 0) {
      Closer(id_);
    }
  }
  HidGuard(const HidGuard&) = delete;
  HidGuard& operator=(const HidGuard&) = delete;

  HidGuard(HidGuard&& other) noexcept : id_(other.id_) { other.id_ = -1; }
  HidGuard& operator=(HidGuard&& other) noexcept {
    id_ = other.id_;
    other.id_ = -1;
    return *this;
  }

  [[nodiscard]] hid_t id() const noexcept { return id_; }
  [[nodiscard]] bool ok() const noexcept { return id_ >= 0; }

  private:
  hid_t id_;
};

using FileGuard = HidGuard<H5Fclose>;
using DsetGuard = HidGuard<H5Dclose>;
using SpaceGuard = HidGuard<H5Sclose>;
using TypeGuard = HidGuard<H5Tclose>;

// ---------- Dimension scale -> coordinate variable name ----------

struct ScaleNameSink {
  std::string name;
  bool found = false;
};

herr_t collectFirstScale(hid_t /*did*/, unsigned /*dim*/, hid_t scaleId, void* op) noexcept {
  auto* sink = static_cast<ScaleNameSink*>(op);
  if (sink->found) {
    return 0;
  }
  char buf[1024]{};
  const ssize_t n = H5Iget_name(scaleId, buf, sizeof(buf));
  if (n > 0) {
    sink->name.assign(buf, static_cast<std::size_t>(n));
    sink->found = true;
    return 1; // stop iteration
  }
  return 0;
}

bool readCoordScalar(hid_t cdset, hsize_t idx, double& out) noexcept {
  const SpaceGuard fspace{H5Dget_space(cdset)};
  if (!fspace.ok()) {
    return false;
  }
  const hsize_t one = 1;
  const SpaceGuard mspace{H5Screate_simple(1, &one, nullptr)};
  if (!mspace.ok()) {
    return false;
  }
  if (H5Sselect_hyperslab(fspace.id(), H5S_SELECT_SET, &idx, nullptr, &one, nullptr) < 0) {
    return false;
  }
  return H5Dread(cdset, H5T_NATIVE_DOUBLE, mspace.id(), fspace.id(), H5P_DEFAULT, &out) >= 0;
}

// ---------- Read distribution ----------

struct DistributeReadCtx {
  hid_t dset;
  void* data;
  std::size_t bytes;
  hid_t memType;
  Grid::Owner owner;
  MPI_Comm comm;
  MPI_Comm nodeComm;
  MPI_Comm leaderComm;
};

/// Performs the H5Dread on whichever rank is supposed to read and, if the
/// broadcast option is active, propagates the result via MPI_Bcast.
bool distributeRead(const DistributeReadCtx& c) noexcept {
  const bool bcast = envOn("ASAGI_LITE_BCAST");

  // Decide who reads and over which communicator the bytes flow.
  MPI_Comm bcastComm = MPI_COMM_NULL;
  bool readsItself = true;

  if (c.owner == Grid::Owner::Shmem) {
    int nodeRank = 0;
    MPI_Comm_rank(c.nodeComm, &nodeRank);
    if (bcast && c.leaderComm != MPI_COMM_NULL) {
      // Only the global leader (rank 0 in leaderComm) reads; result is
      // broadcast to every other node leader. Within-node ranks see it
      // through the shared-memory buffer.
      bcastComm = c.leaderComm;
      int leaderRank = 0;
      MPI_Comm_rank(bcastComm, &leaderRank);
      readsItself = (leaderRank == 0);
    } else if (bcast) {
      // Non-leader rank under broadcast mode: no read, only sync.
      readsItself = false;
    } else {
      // Default shmem path: each node leader reads its own copy.
      readsItself = (nodeRank == 0);
    }
  } else {
    // Owner::Malloc — every rank has its own buffer.
    if (bcast) {
      bcastComm = c.comm;
      int br = 0;
      MPI_Comm_rank(bcastComm, &br);
      readsItself = (br == 0);
    }
    // else: every rank reads independently into its own buffer.
  }

  int readOk = 1;
  if (readsItself) {
    if (H5Dread(c.dset, c.memType, H5S_ALL, H5S_ALL, H5P_DEFAULT, c.data) < 0) {
      readOk = 0;
    }
  }

  if (bcastComm != MPI_COMM_NULL) {
    MPI_Bcast(&readOk, 1, MPI_INT, 0, bcastComm);
    if (readOk != 0) {
      bcastLarge(c.data, c.bytes, 0, bcastComm);
    }
  }

  // Within-node sync: makes the populated shmem buffer visible to every
  // rank on the node (and propagates readOk to non-leaders).
  if (c.owner == Grid::Owner::Shmem) {
    MPI_Bcast(&readOk, 1, MPI_INT, 0, c.nodeComm);
    MPI_Barrier(c.nodeComm);
  }

  return readOk != 0;
}

} // anonymous namespace

// ============================ Lifecycle ============================

Grid::~Grid() noexcept { close(); }

Grid::Grid(Grid&& o) noexcept
    : size_(o.size_), stride_(o.stride_), offset_(o.offset_), invScaling_(o.invScaling_),
      ndims_(o.ndims_), data_(o.data_), elemType_(o.elemType_), owner_(o.owner_), comm_(o.comm_),
      win_(o.win_), nodeComm_(o.nodeComm_), leaderComm_(o.leaderComm_) {
  o.data_ = nullptr;
  o.owner_ = Owner::None;
  o.ndims_ = 0;
  o.win_ = MPI_WIN_NULL;
  o.nodeComm_ = MPI_COMM_NULL;
  o.leaderComm_ = MPI_COMM_NULL;
}

Grid& Grid::operator=(Grid&& o) noexcept {
  if (this != &o) {
    close();
    size_ = o.size_;
    stride_ = o.stride_;
    offset_ = o.offset_;
    invScaling_ = o.invScaling_;
    ndims_ = o.ndims_;
    data_ = o.data_;
    elemType_ = o.elemType_;
    owner_ = o.owner_;
    comm_ = o.comm_;
    win_ = o.win_;
    nodeComm_ = o.nodeComm_;
    leaderComm_ = o.leaderComm_;

    o.data_ = nullptr;
    o.owner_ = Owner::None;
    o.ndims_ = 0;
    o.win_ = MPI_WIN_NULL;
    o.nodeComm_ = MPI_COMM_NULL;
    o.leaderComm_ = MPI_COMM_NULL;
  }
  return *this;
}

std::size_t Grid::sizeAt(std::size_t i) const noexcept { return (i < ndims_) ? size_[i] : 1; }

double Grid::minAt(std::size_t i) const noexcept { return (i < ndims_) ? offset_[i] : 0.0; }

double Grid::maxAt(std::size_t i) const noexcept {
  if (i >= ndims_) {
    return 0.0;
  }
  if (invScaling_[i] == 0.0) {
    return offset_[i];
  }
  return offset_[i] + (static_cast<double>(size_[i]) - 1.0) / invScaling_[i];
}

void Grid::resetMetadata() noexcept {
  ndims_ = 0;
  size_ = {{1, 1, 1, 1}};
  stride_ = {{0, 0, 0, 0}};
  offset_ = {{0.0, 0.0, 0.0, 0.0}};
  invScaling_ = {{1.0, 1.0, 1.0, 1.0}};
}

void Grid::close() noexcept {
  freeBuffer();
  resetMetadata();
}

void Grid::freeBuffer() noexcept {
  switch (owner_) {
  case Owner::Malloc:
    std::free(data_);
    break;
  case Owner::Shmem:
    // All collective on the corresponding communicator.
    if (win_ != MPI_WIN_NULL) {
      MPI_Win_free(&win_);
      win_ = MPI_WIN_NULL;
    }
    if (leaderComm_ != MPI_COMM_NULL) {
      MPI_Comm_free(&leaderComm_);
      leaderComm_ = MPI_COMM_NULL;
    }
    if (nodeComm_ != MPI_COMM_NULL) {
      MPI_Comm_free(&nodeComm_);
      nodeComm_ = MPI_COMM_NULL;
    }
    break;
  case Owner::None:
  default:
    break;
  }
  data_ = nullptr;
  owner_ = Owner::None;
}

// ============================ open() ============================

Error Grid::open(const std::string& filename, const std::string& varname) {
  close();

  const FileGuard fid{H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)};
  if (!fid.ok()) {
    return Error::IoError;
  }

  const DsetGuard dset{H5Dopen2(fid.id(), varname.c_str(), H5P_DEFAULT)};
  if (!dset.ok()) {
    return Error::IoError;
  }

  // -------- Shape --------
  {
    const SpaceGuard space{H5Dget_space(dset.id())};
    const int ndims = H5Sget_simple_extent_ndims(space.id());
    if (ndims < 1 || static_cast<std::size_t>(ndims) > KMaxFileDimensions) {
      return Error::BadRank;
    }
    hsize_t dims[KMaxFileDimensions]{};
    H5Sget_simple_extent_dims(space.id(), dims, nullptr);
    ndims_ = static_cast<std::size_t>(ndims);
    for (std::size_t i = 0; i < ndims_; ++i) {
      size_[i] = static_cast<std::size_t>(dims[i]);
    }
  }

  // -------- Strides (row-major) --------
  stride_[ndims_ - 1] = 1;
  for (std::size_t i = ndims_ - 1; i > 0; --i) {
    stride_[i - 1] = stride_[i] * size_[i];
  }

  // -------- Element type --------
  {
    const TypeGuard dtype{H5Dget_type(dset.id())};
    if (H5Tget_class(dtype.id()) != H5T_FLOAT) {
      return Error::BadElementType;
    }
    const std::size_t tsize = H5Tget_size(dtype.id());
    if (tsize == 4) {
      elemType_ = ElementType::Float;
    } else if (tsize == 8) {
      elemType_ = ElementType::Double;
    } else {
      return Error::BadElementType;
    }
  }

  // -------- Tier decision + buffer allocation --------
  std::size_t nelem = 1;
  for (std::size_t i = 0; i < ndims_; ++i) {
    nelem *= size_[i];
  }
  const std::size_t elemSize = (elemType_ == ElementType::Float) ? 4 : 8;
  const std::size_t bytes = nelem * elemSize;

  MPI_Comm probeNodeComm = MPI_COMM_NULL;
  const bool useShmem = decideShmem(bytes, comm_, probeNodeComm);

  if (useShmem) {
    // Hand the already-split node communicator over to the member.
    nodeComm_ = probeNodeComm;
    if (auto e = allocateShared(bytes); e != Error::Success) {
      close();
      return e;
    }
  } else {
    // Probe communicator is unused; free it.
    if (probeNodeComm != MPI_COMM_NULL) {
      MPI_Comm_free(&probeNodeComm);
    }
    if (auto e = allocateLocal(bytes); e != Error::Success) {
      close();
      return e;
    }
  }

  // -------- Read + (optional) broadcast --------
  {
    const DistributeReadCtx ctx{
        /*dset=*/dset.id(),
        /*data=*/data_,
        /*bytes=*/bytes,
        /*memType=*/(elemType_ == ElementType::Float) ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE,
        /*owner=*/owner_,
        /*comm=*/comm_,
        /*nodeComm=*/nodeComm_,
        /*leaderComm=*/leaderComm_,
    };
    if (!distributeRead(ctx)) {
      close();
      return Error::IoError;
    }
  }

  // -------- Coordinate variables -> offset/scaling --------
  for (std::size_t i = 0; i < ndims_; ++i) {
    ScaleNameSink sink;
    H5DSiterate_scales(dset.id(), static_cast<unsigned>(i), nullptr, collectFirstScale, &sink);

    if (!sink.found) {
      offset_[i] = 0.0;
      invScaling_[i] = 1.0;
      continue;
    }
    const DsetGuard coord{H5Dopen2(fid.id(), sink.name.c_str(), H5P_DEFAULT)};
    if (!coord.ok()) {
      offset_[i] = 0.0;
      invScaling_[i] = 1.0;
      continue;
    }
    double v0 = 0.0;
    double vN = 0.0;
    if (!readCoordScalar(coord.id(), 0, v0)) {
      offset_[i] = 0.0;
      invScaling_[i] = 1.0;
      continue;
    }
    if (size_[i] <= 1) {
      offset_[i] = v0;
      invScaling_[i] = 0.0;
    } else if (readCoordScalar(coord.id(), size_[i] - 1, vN) && vN != v0) {
      offset_[i] = v0;
      invScaling_[i] = static_cast<double>(size_[i] - 1) / (vN - v0);
    } else {
      offset_[i] = v0;
      invScaling_[i] = 0.0;
    }
  }

  return Error::Success;
}

// ============================ Pointer variants ============================

void Grid::coordsToIndices(const double* coords,
                           const Rounding* rounding,
                           std::size_t n,
                           std::size_t* outIndices) const noexcept {
  assert(n <= ndims_ && "more coordinates than grid dimensions");
  for (std::size_t i = 0; i < n; ++i) {
    const double pos = coordPos(i, coords[i]);
    const Rounding r = (rounding != nullptr) ? rounding[i] : Rounding::Nearest;
    switch (r) {
    case Rounding::Down:
      outIndices[i] = static_cast<std::size_t>(std::floor(pos));
      break;
    case Rounding::Up:
      outIndices[i] = static_cast<std::size_t>(std::ceil(pos));
      break;
    case Rounding::Nearest:
    default:
      outIndices[i] = static_cast<std::size_t>(std::lround(pos));
      break;
    }
  }
}

template <typename T>
void Grid::getValue(const std::size_t* indices,
                    std::size_t nIndices,
                    T* out,
                    std::size_t mComponents) const noexcept {
  static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "T must be float or double");

  assert(data_ != nullptr && "grid is not open");
  assert(((mComponents == 1 && ndims_ == nIndices) ||
          (ndims_ == nIndices + 1 && size_[nIndices] == mComponents)) &&
         "grid shape does not match (nIndices, mComponents)");

  std::size_t idx = 0;
  for (std::size_t i = 0; i < nIndices; ++i) {
    assert(indices[i] < size_[i] && "index out of range");
    idx += indices[i] * stride_[i];
  }

  if (elemType_ == ElementType::Float) {
    const float* p = static_cast<const float*>(data_) + idx;
    for (std::size_t j = 0; j < mComponents; ++j) {
      out[j] = static_cast<T>(p[j]);
    }
  } else {
    const double* p = static_cast<const double*>(data_) + idx;
    for (std::size_t j = 0; j < mComponents; ++j) {
      out[j] = static_cast<T>(p[j]);
    }
  }
}

template void
    Grid::getValue<float>(const std::size_t*, std::size_t, float*, std::size_t) const noexcept;
template void
    Grid::getValue<double>(const std::size_t*, std::size_t, double*, std::size_t) const noexcept;

// ============================ Allocation ============================

Error Grid::allocateLocal(std::size_t bytes) {
  void* p = std::malloc(bytes);
  if (p == nullptr) {
    return Error::AllocationFailed;
  }
  data_ = p;
  owner_ = Owner::Malloc;
  return Error::Success;
}

Error Grid::allocateShared(std::size_t bytes) {
  // nodeComm_ has already been set up by the caller (decideShmem()).
  int nodeRank = 0;
  MPI_Comm_rank(nodeComm_, &nodeRank);

  const MPI_Aint ourBytes = (nodeRank == 0) ? static_cast<MPI_Aint>(bytes) : MPI_Aint{0};
  void* base = nullptr;
  const int dispUnit = 1; // bytes
  if (MPI_Win_allocate_shared(
          ourBytes, dispUnit, MPI_INFO_NULL, nodeComm_, static_cast<void*>(&base), &win_) !=
      MPI_SUCCESS) {
    return Error::AllocationFailed;
  }
  if (nodeRank != 0) {
    MPI_Aint qSize = 0;
    int qDisp = 0;
    MPI_Win_shared_query(win_, 0, &qSize, &qDisp, static_cast<void*>(&base));
  }
  data_ = base;
  owner_ = Owner::Shmem;

  // Build a leader-only sub-communicator (one rank per node) when broadcast
  // is requested. Required so the global reader can address peer leaders
  // directly without involving the within-node ranks.
  if (envOn("ASAGI_LITE_BCAST")) {
    int globalRank = 0;
    MPI_Comm_rank(comm_, &globalRank);
    const int color = (nodeRank == 0) ? 0 : MPI_UNDEFINED;
    if (MPI_Comm_split(comm_, color, globalRank, &leaderComm_) != MPI_SUCCESS) {
      // Non-fatal: degrade to per-leader-reads behaviour.
      leaderComm_ = MPI_COMM_NULL;
    }
  }
  return Error::Success;
}

} // namespace asagi_lite
