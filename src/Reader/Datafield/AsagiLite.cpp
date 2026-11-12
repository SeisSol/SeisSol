// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "AsagiLite.h"

#include "Grid.h"
#include "Interpolation.h"
#include "utils/logger.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <mpi.h>
#include <new>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sys/types.h>
#include <unistd.h> // sysconf
#include <utility>
#include <vector>

namespace seissol::reader::datafield {
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
  return bytes >= AsagiLiteGrid::SharedMemThreshold;
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

/// The string attribute `name` of `object`, of fixed or variable length; empty where it has none.
std::string stringAttribute(hid_t object, const char* name) {
  if (H5Aexists(object, name) <= 0) {
    return {};
  }
  const hid_t attribute = H5Aopen(object, name, H5P_DEFAULT);
  if (attribute < 0) {
    return {};
  }
  std::string value;
  const TypeGuard type{H5Aget_type(attribute)};
  if (type.ok() && H5Tget_class(type.id()) == H5T_STRING) {
    // in the character set of the file, since HDF5 converts between none
    const TypeGuard memory{H5Tcopy(H5T_C_S1)};
    H5Tset_cset(memory.id(), H5Tget_cset(type.id()));
    if (H5Tis_variable_str(type.id()) > 0) {
      char* data = nullptr;
      H5Tset_size(memory.id(), H5T_VARIABLE);
      if (H5Aread(attribute, memory.id(), static_cast<void*>(&data)) >= 0 && data != nullptr) {
        value = data;
        H5free_memory(data);
      }
    } else {
      const std::size_t size = H5Tget_size(type.id());
      std::string buffer(size + 1, '\0');
      H5Tset_size(memory.id(), size + 1);
      if (H5Aread(attribute, memory.id(), buffer.data()) >= 0) {
        value = buffer.c_str();
      }
    }
  }
  H5Aclose(attribute);
  return value;
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

/// The whole coordinate variable `cdset` of `count` values, or nothing where it has another
/// number of values or cannot be read.
std::optional<std::vector<double>> readCoordinates(hid_t cdset, std::size_t count) noexcept {
  const SpaceGuard fspace{H5Dget_space(cdset)};
  if (!fspace.ok() || H5Sget_simple_extent_ndims(fspace.id()) != 1) {
    return std::nullopt;
  }
  hsize_t extent = 0;
  H5Sget_simple_extent_dims(fspace.id(), &extent, nullptr);
  if (static_cast<std::size_t>(extent) != count) {
    return std::nullopt;
  }
  try {
    std::vector<double> values(count);
    if (H5Dread(cdset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0) {
      return std::nullopt;
    }
    return values;
  } catch (const std::bad_alloc&) {
    return std::nullopt;
  }
}

// ---------- Read distribution ----------

struct DistributeReadCtx {
  hid_t dset;
  void* data;
  std::size_t bytes;
  hid_t memType;
  /// H5S_ALL for a whole-dataset read, or a selected file/memory space pair for
  /// a window slice. `data` and `bytes` describe only the region being read, so
  /// the broadcast moves the slice and not the whole buffer.
  hid_t fileSpace;
  hid_t memSpace;
  AsagiLiteGrid::Owner owner;
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

  if (c.owner == AsagiLiteGrid::Owner::Shmem) {
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
    if (H5Dread(c.dset, c.memType, c.memSpace, c.fileSpace, H5P_DEFAULT, c.data) < 0) {
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
  if (c.owner == AsagiLiteGrid::Owner::Shmem) {
    MPI_Bcast(&readOk, 1, MPI_INT, 0, c.nodeComm);
    MPI_Barrier(c.nodeComm);
  }

  return readOk != 0;
}

} // anonymous namespace

std::optional<std::vector<AsagiAxis>> asagiAxesOf(const std::string& file,
                                                  const std::string& variable) {
  const FileGuard handle{H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)};
  if (!handle.ok()) {
    return std::nullopt;
  }
  const DsetGuard dataset{H5Dopen2(handle.id(), variable.c_str(), H5P_DEFAULT)};
  if (!dataset.ok()) {
    return std::nullopt;
  }
  const SpaceGuard space{H5Dget_space(dataset.id())};
  const int ndims = H5Sget_simple_extent_ndims(space.id());
  if (ndims < 1 || static_cast<std::size_t>(ndims) > AsagiLiteGrid::MaxFileDimensions) {
    return std::nullopt;
  }
  hsize_t dims[AsagiLiteGrid::MaxFileDimensions]{};
  H5Sget_simple_extent_dims(space.id(), dims, nullptr);
  const TypeGuard type{H5Dget_type(dataset.id())};
  const bool compound = H5Tget_class(type.id()) == H5T_COMPOUND;
  std::vector<AsagiAxis> axes;
  for (int i = 0; i < ndims; ++i) {
    ScaleNameSink sink;
    H5DSiterate_scales(dataset.id(), static_cast<unsigned>(i), nullptr, collectFirstScale, &sink);
    const auto nodes = static_cast<std::size_t>(dims[i]);
    // the component axis of a vector grid, as open() tells it
    if (!sink.found && i + 1 == ndims && ndims > 1 && !compound) {
      continue;
    }
    AsagiAxis axis{0.0, static_cast<double>(nodes) - 1.0, nodes};
    if (sink.found) {
      const DsetGuard coordinate{H5Dopen2(handle.id(), sink.name.c_str(), H5P_DEFAULT)};
      if (coordinate.ok()) {
        readCoordScalar(coordinate.id(), 0, axis.first);
        readCoordScalar(coordinate.id(), nodes - 1, axis.last);
      }
    }
    axes.push_back(axis);
  }
  std::reverse(axes.begin(), axes.end());
  return axes;
}

// ============================ Lifecycle ============================

AsagiLiteGrid::~AsagiLiteGrid() noexcept { close(); }

// Move operations are gone: a Grid is owned through unique_ptr by GridStore and
// never relocated, and a moved-from grid holding a live MPI window is a hazard
// with no upside.

void AsagiLiteGrid::resetMetadata() noexcept {
  ndims_ = 0;
  size_ = {{1, 1, 1, 1}};
  stride_ = {{0, 0, 0, 0}};
  offset_ = {{0.0, 0.0, 0.0, 0.0}};
  deltaInv_ = {{1.0, 1.0, 1.0, 1.0}};
  delta_ = {{1.0, 1.0, 1.0, 1.0}};
  hasComponentAxis_ = false;
  members_ = 0;
  window_ = TimeWindow{};
  numTimeSlices_ = 1;
  sliceElems_ = 0;
  sliceBytes_ = 0;
  timeOrigin_ = 0.0;
  for (auto& knots : knots_) {
    knots.clear();
  }
  for (auto& spacings : inverseSpacings_) {
    spacings.clear();
  }
  view_ = ArrayView{};
  ++revision_;
}

void AsagiLiteGrid::close() noexcept {
  freeBuffer();
  closeFile();
  resetMetadata();
}

void AsagiLiteGrid::closeFile() noexcept {
  // a grid of a store that lives until the process ends outlives HDF5, which closes every handle
  // when it terminates: then there is nothing left to close
  const auto close = [](hid_t& handle, herr_t (*closer)(hid_t)) {
    if (handle >= 0 && H5Iis_valid(handle) > 0) {
      closer(handle);
    }
    handle = -1;
  };
  close(memType_, H5Tclose);
  close(dsetId_, H5Dclose);
  close(fileId_, H5Fclose);
}

void AsagiLiteGrid::freeBuffer() noexcept {
  int finalized = 0;
  MPI_Finalized(&finalized);
  switch (owner_) {
  case Owner::Malloc:
    std::free(data_);
    break;
  case Owner::Shmem:
    // All collective on the corresponding communicator; after MPI_Finalize -- a grid of a store
    // that lives until the process ends -- the end of the process frees them.
    if (finalized != 0) {
      win_ = MPI_WIN_NULL;
      leaderComm_ = MPI_COMM_NULL;
      nodeComm_ = MPI_COMM_NULL;
      break;
    }
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

AsagiLiteGrid::Error AsagiLiteGrid::open(const std::string& filename,
                                         const std::string& varname,
                                         Interpolation interp,
                                         BoundaryMode boundary,
                                         std::optional<int> timeAxis,
                                         const std::vector<AxisKind>& axes) {
  interp_ = interp;
  boundary_ = boundary;
  timeAxis_ = timeAxis;
  path_ = filename;
  close();
  axes_ = axes;

  // The handles stay open: a sliding window reads again on every advance.
  fileId_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (fileId_ < 0) {
    return Error::IoError;
  }
  dsetId_ = H5Dopen2(fileId_, varname.c_str(), H5P_DEFAULT);
  if (dsetId_ < 0) {
    close();
    return Error::IoError;
  }
  const hid_t fidRaw = fileId_;
  const hid_t dsetRaw = dsetId_;

  // -------- Shape --------
  {
    const SpaceGuard space{H5Dget_space(dsetRaw)};
    const int ndims = H5Sget_simple_extent_ndims(space.id());
    if (ndims < 1 || static_cast<std::size_t>(ndims) > MaxFileDimensions) {
      return Error::BadRank;
    }
    hsize_t dims[MaxFileDimensions]{};
    H5Sget_simple_extent_dims(space.id(), dims, nullptr);
    ndims_ = static_cast<std::size_t>(ndims);
    for (std::size_t i = 0; i < ndims_; ++i) {
      size_[i] = static_cast<std::size_t>(dims[i]);
    }
  }

  // -------- Element type --------
  //
  // Floating point of 32 or 64 bits; or a compound of such members, all of one size -- the form
  // of ASAGI's files (`compound material { float rho; float mu; } data(y, x)`) -- whose members
  // are the components, in their order, read side by side as a trailing component axis would be.
  {
    const TypeGuard dtype{H5Dget_type(dsetRaw)};
    std::size_t tsize = 0;
    if (H5Tget_class(dtype.id()) == H5T_COMPOUND) {
      const int members = H5Tget_nmembers(dtype.id());
      if (members < 1) {
        return Error::BadElementType;
      }
      for (int m = 0; m < members; ++m) {
        const TypeGuard member{H5Tget_member_type(dtype.id(), static_cast<unsigned>(m))};
        const std::size_t size = H5Tget_size(member.id());
        if (H5Tget_class(member.id()) != H5T_FLOAT || (tsize != 0 && size != tsize)) {
          return Error::BadElementType;
        }
        tsize = size;
      }
      members_ = static_cast<std::size_t>(members);
    } else if (H5Tget_class(dtype.id()) == H5T_FLOAT) {
      tsize = H5Tget_size(dtype.id());
    } else {
      return Error::BadElementType;
    }
    if (tsize == 4) {
      elemType_ = ElementType::Float;
    } else if (tsize == 8) {
      elemType_ = ElementType::Double;
    } else {
      return Error::BadElementType;
    }
    if (members_ > 0) {
      // the members side by side, of the native type, matched by name
      memType_ = H5Tcreate(H5T_COMPOUND, members_ * tsize);
      for (std::size_t m = 0; m < members_; ++m) {
        char* member = H5Tget_member_name(dtype.id(), static_cast<unsigned>(m));
        H5Tinsert(memType_,
                  member,
                  m * tsize,
                  elemType_ == ElementType::Float ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE);
        H5free_memory(member);
      }
    }
  }

  // -------- Strides (row-major), in values: the members of a compound side by side --------
  stride_[ndims_ - 1] = std::max<std::size_t>(members_, 1);
  for (std::size_t i = ndims_ - 1; i > 0; --i) {
    stride_[i - 1] = stride_[i] * size_[i];
  }

  // -------- Coordinate variables -> offset/delta --------
  //
  // `delta` is the mean spacing between the ends, (v[n-1] - v[0]) / (n - 1), and the
  // reciprocal is formed as 1.0 / delta: ASAGI's derivation (NetCdfReader::getScaling()), with
  // which easi samples the file, and matching it is the whole point. See the geometry-derivation
  // note in Grid.h.
  fileAxes_.assign(ndims_, GridAxisInfo{});
  // per file dimension its coordinate variable, and the extent its first spacing implies where
  // that is not the extent of the file
  std::vector<std::string> scales(ndims_);
  std::vector<std::optional<long>> nonUniform(ndims_);
  for (std::size_t i = 0; i < ndims_; ++i) {
    ScaleNameSink sink;
    H5DSiterate_scales(dsetRaw, static_cast<unsigned>(i), nullptr, collectFirstScale, &sink);
    if (sink.found) {
      fileAxes_[i].dimension =
          sink.name.size() > 1 && sink.name[0] == '/' ? sink.name.substr(1) : sink.name;
      scales[i] = sink.name;
    }

    // A trailing dimension with no dimension scale is the component axis of a
    // vector grid -- there is no coordinate to interpolate along it. This is the
    // only signal in the file that tells [d_0..d_{N-1}, M] apart from an
    // (N+1)-dimensional scalar grid.
    if (!sink.found && i + 1 == ndims_ && ndims_ > 1 && members_ == 0) {
      hasComponentAxis_ = true;
      offset_[i] = 0.0;
      delta_[i] = 1.0;
      deltaInv_[i] = 1.0;
      continue;
    }

    const auto degenerate = [&](double origin) {
      offset_[i] = origin;
      delta_[i] = 0.0;
      deltaInv_[i] = 0.0;
    };
    const auto unitless = [&]() {
      offset_[i] = 0.0;
      delta_[i] = 1.0;
      deltaInv_[i] = 1.0;
    };

    if (!sink.found) {
      unitless();
      continue;
    }
    const DsetGuard coord{H5Dopen2(fidRaw, sink.name.c_str(), H5P_DEFAULT)};
    if (!coord.ok()) {
      unitless();
      continue;
    }
    fileAxes_[i].units = stringAttribute(coord.id(), "units");
    fileAxes_[i].positive = stringAttribute(coord.id(), "positive");
    fileAxes_[i].standardName = stringAttribute(coord.id(), "standard_name");
    double v0 = 0.0;
    double v1 = 0.0;
    if (!readCoordScalar(coord.id(), 0, v0)) {
      unitless();
      continue;
    }
    if (size_[i] <= 1) {
      degenerate(v0);
      continue;
    }
    double vN = 0.0;
    if (!readCoordScalar(coord.id(), size_[i] - 1, vN) || !(std::isfinite(vN - v0)) || vN == v0) {
      degenerate(v0);
      continue;
    }
    // ASAGI's scaling, with which easi samples the file: the mean spacing between the ends
    offset_[i] = v0;
    delta_[i] = (vN - v0) / static_cast<double>(size_[i] - 1);
    deltaInv_[i] = 1.0 / delta_[i];

    // Cross-check the uniformity assumption: if the first spacing does not reach the far end
    // within half a spacing, the axis is not uniform and the uniform-grid model does not hold for
    // this file. Keep the HDF5 extent -- it is known exactly -- and say so, rather than silently
    // interpolating on a geometry that is not there.
    if (readCoordScalar(coord.id(), 1, v1) && std::isfinite(v1 - v0) && v1 != v0) {
      const long reconstructed = std::lround((vN - v0) / (v1 - v0)) + 1;
      if (reconstructed != static_cast<long>(size_[i])) {
        nonUniform[i] = reconstructed;
      }
    }
  }

  // -------- Axes of knots --------
  //
  // An axis declared knots reads its coordinate variable whole: a lookup searches it. The others
  // stay uniform, from the mean spacing between their ends; one that is not is worth a warning.
  const std::size_t nspatial = ndims_ - (hasComponentAxis_ ? 1 : 0);
  if (!axes_.empty() && axes_.size() != nspatial) {
    logWarning() << "AsagiLite:" << filename << "declares" << axes_.size() << "axes of" << varname
                 << ", which has" << nspatial << ".";
    close();
    return Error::BadKnots;
  }
  for (std::size_t d = 0; d < nspatial; ++d) {
    const std::size_t f = nspatial - 1 - d;
    if (d >= axes_.size() || axes_[d] != AxisKind::Knots) {
      if (nonUniform[f].has_value()) {
        logWarning() << "AsagiLite:" << filename << "axis" << f << "of" << varname
                     << "is not uniform: extent is" << size_[f]
                     << "points but the first spacing implies" << *nonUniform[f]
                     << ". Sampling will use the mean spacing between its ends and the true"
                     << "extent; an axis of knots (axes in its declaration) samples it as it is.";
      }
      continue;
    }
    const auto problem = [&](const std::string& what) {
      logWarning() << "AsagiLite:" << filename << "declares axis" << d << "of" << varname
                   << "knots, but" << what << ".";
      close();
      return Error::BadKnots;
    };
    if (timeAxis_.has_value() && static_cast<std::size_t>(*timeAxis_) == d) {
      return problem("it is the time axis, which is uniform");
    }
    if (size_[f] > MaxKnotsOfFile) {
      return problem("it has " + std::to_string(size_[f]) + " knots, more than " +
                     std::to_string(MaxKnotsOfFile));
    }
    if (scales[f].empty()) {
      return problem("the file has no coordinate variable for it");
    }
    const DsetGuard coord{H5Dopen2(fidRaw, scales[f].c_str(), H5P_DEFAULT)};
    std::optional<std::vector<double>> knots;
    if (coord.ok()) {
      knots = readCoordinates(coord.id(), size_[f]);
    }
    if (!knots.has_value()) {
      return problem("its coordinate variable " + scales[f] + " cannot be read as " +
                     std::to_string(size_[f]) + " values");
    }
    for (std::size_t k = 0; k < knots->size(); ++k) {
      if (!std::isfinite((*knots)[k])) {
        return problem("its knot " + std::to_string(k) + " is not finite");
      }
      if (k > 0 && !((*knots)[k - 1] < (*knots)[k])) {
        return problem("its knots do not increase strictly, at knot " + std::to_string(k));
      }
    }
    inverseSpacings_[d].clear();
    for (std::size_t k = 0; k + 1 < knots->size(); ++k) {
      inverseSpacings_[d].push_back(1.0 / ((*knots)[k + 1] - (*knots)[k]));
    }
    // a lookup reads the spacing of its cell, the first of a single knot as well
    if (inverseSpacings_[d].empty()) {
      inverseSpacings_[d].push_back(0.0);
    }
    knots_[d] = std::move(*knots);
  }

  buildView();

  // -------- Time axis --------
  //
  // A window is only a contiguous hyperslab if time is the slowest file
  // dimension. Anything else would turn every slide into a strided gather, so
  // it is rejected rather than silently paid for.
  const std::size_t elemSize = (elemType_ == ElementType::Float) ? 4 : 8;
  std::size_t nelem = std::max<std::size_t>(members_, 1);
  for (std::size_t i = 0; i < ndims_; ++i) {
    nelem *= size_[i];
  }

  if (timeAxis_.has_value()) {
    const int axis = *timeAxis_;
    if (axis < 0 || static_cast<std::size_t>(axis) + 1 != nspatial) {
      logError() << "AsagiLite:" << filename << "declares grid axis" << axis
                 << "as its time axis, but only axis" << (static_cast<long>(nspatial) - 1)
                 << "can be — the time axis must be the slowest dimension in the file"
                 << "(the NetCDF convention, z(t, y, x)) so that a resident window is one"
                 << "contiguous hyperslab.";
      close();
      return Error::BadRank;
    }
    numTimeSlices_ = size_[0];
    sliceElems_ = nelem / size_[0];
    timeOrigin_ = offset_[0];
  } else {
    numTimeSlices_ = 1;
    sliceElems_ = nelem;
    timeOrigin_ = 0.0;
  }
  sliceBytes_ = sliceElems_ * elemSize;

  // A static grid has nothing to size, so load it now and be done.
  if (!timeAxis_.has_value()) {
    return resizeWindow(1);
  }
  return Error::Success;
}

// ============================ Window management ============================

AsagiLiteGrid::Error AsagiLiteGrid::resizeWindow(std::size_t residentSlices) {
  const std::size_t count = std::min(std::max(residentSlices, std::size_t{1}), numTimeSlices_);
  const std::size_t bytes = count * sliceBytes_;

  freeBuffer();
  MPI_Comm probeNodeComm = MPI_COMM_NULL;
  const bool useShmem = decideShmem(bytes, comm_, probeNodeComm);
  if (useShmem) {
    nodeComm_ = probeNodeComm;
    if (auto e = allocateShared(bytes); e != Error::Success) {
      close();
      return e;
    }
  } else {
    if (probeNodeComm != MPI_COMM_NULL) {
      MPI_Comm_free(&probeNodeComm);
    }
    if (auto e = allocateLocal(bytes); e != Error::Success) {
      close();
      return e;
    }
  }

  window_ = TimeWindow{0, count};
  if (auto e = readSlices(0, count, 0); e != Error::Success) {
    close();
    return e;
  }
  buildView();
  return Error::Success;
}

AsagiLiteGrid::Error AsagiLiteGrid::advanceTo(double time) {
  if (!timeAxis_.has_value() || window_.count >= numTimeSlices_) {
    return Error::Success; // fully resident: nothing can expire
  }
  const auto axis = static_cast<unsigned>(*timeAxis_);
  const double deltaInv = deltaInv_[0];
  const double raw = (time - timeOrigin_) * deltaInv;
  std::int64_t base = 0;
  if (raw > 0.0) {
    base = static_cast<std::int64_t>(std::floor(raw));
    const auto maxBase = static_cast<std::int64_t>(numTimeSlices_) - 2;
    if (base > maxBase) {
      base = maxBase > 0 ? maxBase : 0;
    }
  }

  // Placed anew whenever `time` moved on by a cell, not only where the window no longer serves
  // `time` itself: the window has to serve every time until the next update, which comes up to
  // GridStore::suggestedSyncInterval() later, (W - w - 2s) cells -- exactly the span ahead of a
  // window placed so (placeTimeWindow()), with s = TimeWindowSlack slices to spare on either side.
  // A window kept while it still serves `time` may serve only a cell more, and the samples after
  // it would be clamped onto its last slice.
  const StencilKernel kernel = kernelOf(interp_);
  const TimeWindow target = placeTimeWindow(base, kernel, window_.count, numTimeSlices_);
  if (target == window_) {
    return Error::Success;
  }

  // Reuse the overlap. Only the buffer owner may touch it: with shared memory
  // that is the node leader, and the barriers keep the other ranks out of a
  // buffer that is mid-move.
  const std::int64_t shift =
      static_cast<std::int64_t>(target.start) - static_cast<std::int64_t>(window_.start);
  const auto count = static_cast<std::int64_t>(window_.count);
  const std::int64_t overlap = (std::abs(shift) < count) ? count - std::abs(shift) : 0;

  bool owns = true;
  if (owner_ == Owner::Shmem) {
    int nodeRank = 0;
    MPI_Comm_rank(nodeComm_, &nodeRank);
    owns = (nodeRank == 0);
    MPI_Barrier(nodeComm_);
  }
  if (overlap > 0 && owns) {
    auto* bytes = static_cast<char*>(data_);
    const std::size_t src = (shift > 0) ? sliceOffset(static_cast<std::size_t>(shift)) : 0;
    const std::size_t dst = (shift > 0) ? 0 : sliceOffset(static_cast<std::size_t>(-shift));
    const std::size_t elemSize = (elemType_ == ElementType::Float) ? 4 : 8;
    std::memmove(bytes + dst * elemSize,
                 bytes + src * elemSize,
                 static_cast<std::size_t>(overlap) * sliceBytes_);
  }
  if (owner_ == Owner::Shmem) {
    MPI_Barrier(nodeComm_);
  }

  window_ = target;
  Error err = Error::Success;
  if (overlap == 0) {
    err = readSlices(target.start, target.count, 0);
  } else if (shift > 0) {
    err = readSlices(target.start + static_cast<std::size_t>(overlap),
                     static_cast<std::size_t>(shift),
                     static_cast<std::size_t>(overlap));
  } else {
    err = readSlices(target.start, static_cast<std::size_t>(-shift), 0);
  }
  if (err != Error::Success) {
    return err;
  }
  // view_.min on the time axis tracks the window, not the file.
  view_.min[axis] = timeOrigin_ + static_cast<double>(target.start) * delta_[0];
  view_.before[axis] = target.start;
  view_.after[axis] = numTimeSlices_ - target.start - target.count;
  ++revision_;
  return Error::Success;
}

AsagiLiteGrid::Error
    AsagiLiteGrid::readSlices(std::size_t first, std::size_t count, std::size_t destSlot) {
  if (count == 0) {
    return Error::Success;
  }
  const std::size_t elemSize = (elemType_ == ElementType::Float) ? 4 : 8;
  void* dest = static_cast<char*>(data_) + sliceOffset(destSlot) * elemSize;
  const hid_t memType = members_ > 0                        ? memType_
                        : (elemType_ == ElementType::Float) ? H5T_NATIVE_FLOAT
                                                            : H5T_NATIVE_DOUBLE;

  hid_t fileSpace = H5S_ALL;
  hid_t memSpace = H5S_ALL;
  SpaceGuard fspace{-1};
  SpaceGuard mspace{-1};
  const bool partial = timeAxis_.has_value() && count < numTimeSlices_;
  if (partial) {
    hsize_t start[MaxFileDimensions]{};
    hsize_t block[MaxFileDimensions]{};
    start[0] = first;
    block[0] = count;
    for (std::size_t i = 1; i < ndims_; ++i) {
      start[i] = 0;
      block[i] = size_[i];
    }
    fspace = SpaceGuard{H5Dget_space(dsetId_)};
    if (!fspace.ok() ||
        H5Sselect_hyperslab(fspace.id(), H5S_SELECT_SET, start, nullptr, block, nullptr) < 0) {
      return Error::IoError;
    }
    mspace = SpaceGuard{H5Screate_simple(static_cast<int>(ndims_), block, nullptr)};
    if (!mspace.ok()) {
      return Error::IoError;
    }
    fileSpace = fspace.id();
    memSpace = mspace.id();
  }

  const DistributeReadCtx ctx{
      /*dset=*/dsetId_,
      /*data=*/dest,
      /*bytes=*/count * sliceBytes_,
      /*memType=*/memType,
      /*fileSpace=*/fileSpace,
      /*memSpace=*/memSpace,
      /*owner=*/owner_,
      /*comm=*/comm_,
      /*nodeComm=*/nodeComm_,
      /*leaderComm=*/leaderComm_,
  };
  return distributeRead(ctx) ? Error::Success : Error::IoError;
}

// ============================ Grid interface ============================

// Permutes the file-order metadata into grid order, once. Grid axis d is file
// dimension (nspatial - 1 - d), with the component axis stripped first so it
// never appears as a spatial axis.
void AsagiLiteGrid::buildView() noexcept {
  const std::size_t nspatial = ndims_ - (hasComponentAxis_ ? 1 : 0);
  ++revision_;

  view_ = ArrayView{};
  view_.data = data_;
  view_.type = elemType_;
  view_.dimensions = static_cast<unsigned>(nspatial);
  view_.components = members_ > 0        ? static_cast<unsigned>(members_)
                     : hasComponentAxis_ ? static_cast<unsigned>(size_[ndims_ - 1])
                                         : 1U;
  view_.componentStride = hasComponentAxis_ ? stride_[ndims_ - 1] : 1;

  for (std::size_t d = 0; d < nspatial; ++d) {
    const std::size_t f = nspatial - 1 - d;
    view_.min[d] = offset_[f];
    view_.deltaInv[d] = deltaInv_[f];
    view_.num[d] = static_cast<unsigned>(size_[f]);
    view_.stride[d] = stride_[f];
    if (d < knots_.size() && !knots_[d].empty()) {
      view_.min[d] = knots_[d].front();
      view_.knots[d] = knots_[d].data();
      view_.inverseSpacings[d] = inverseSpacings_[d].data();
    }
  }
  // On the time axis the view sees only the resident window, and its origin is
  // the window's first slice rather than the file's. Everything downstream —
  // sampling, bounds(), the stencil edge fallback — then works on the window
  // without knowing a window exists.
  if (timeAxis_.has_value() && window_.count > 0) {
    const auto axis = static_cast<unsigned>(*timeAxis_);
    if (axis < nspatial) {
      view_.num[axis] = static_cast<unsigned>(window_.count);
      view_.min[axis] = timeOrigin_ + static_cast<double>(window_.start) * delta_[0];
      view_.before[axis] = window_.start;
      view_.after[axis] = numTimeSlices_ - window_.start - window_.count;
    }
  }
}

std::vector<GridAxisInfo> AsagiLiteGrid::axisInfo() const {
  // the grid axes are the spatial dimensions of the file, in the opposite order
  const std::size_t nspatial = view_.dimensions;
  std::vector<GridAxisInfo> axes(nspatial);
  for (std::size_t d = 0; d < nspatial && nspatial - 1 - d < fileAxes_.size(); ++d) {
    axes[d] = fileAxes_[nspatial - 1 - d];
  }
  return axes;
}

void AsagiLiteGrid::geometry(GridGeometry& out) const {
  out = GridGeometry{};
  out.dimensions = view_.dimensions;
  const std::size_t nspatial = view_.dimensions;
  for (std::size_t d = 0; d < nspatial; ++d) {
    out.min[d] = view_.min[d];
    out.delta[d] = delta_[nspatial - 1 - d];
    out.num[d] = view_.num[d];
    // knots have no one spacing: the mean, which spans the axis as the knots do
    if (d < knots_.size() && knots_[d].size() > 1) {
      out.delta[d] =
          (knots_[d].back() - knots_[d].front()) / static_cast<double>(knots_[d].size() - 1);
    }
  }
}

bool AsagiLiteGrid::lookupsOnly() const noexcept {
  return lookupsOnly_ || std::any_of(knots_.begin(), knots_.end(), [](const auto& knots) {
           return !knots.empty();
         });
}

void AsagiLiteGrid::sample(const double* x, double* out) const {
  if (lookupsOnly()) {
    throw std::invalid_argument(
        "datafield: a grid with axes of knots, or an interpolation or a "
        "boundary per axis, is sampled through the written-out lookup only");
  }
  const std::size_t clamped = datafield::sample(view_, interp_, x, out);
  if (clamped != 0 && boundary_ == BoundaryMode::Fail) {
    logError() << "AsagiLite:" << path_ << "queried outside its volume.";
  }
}

void AsagiLiteGrid::sampleBatch(const double* const* x,
                                std::size_t count,
                                double* const* out) const {
  if (lookupsOnly()) {
    throw std::invalid_argument(
        "datafield: a grid with axes of knots, or an interpolation or a "
        "boundary per axis, is sampled through the written-out lookup only");
  }
  const std::size_t clamped = datafield::sampleBatch(view_, interp_, x, count, out, scratch_);
  reportClamped(clamped, count);
}

void AsagiLiteGrid::sampleBatch(const double* const* x,
                                std::size_t count,
                                float* const* out) const {
  if (lookupsOnly()) {
    throw std::invalid_argument(
        "datafield: a grid with axes of knots, or an interpolation or a "
        "boundary per axis, is sampled through the written-out lookup only");
  }
  const std::size_t clamped = datafield::sampleBatch(view_, interp_, x, count, out, scratch_);
  reportClamped(clamped, count);
}

void AsagiLiteGrid::reportClamped(std::size_t clamped, std::size_t count) const {
  if (clamped != 0 && boundary_ == BoundaryMode::Fail) {
    logError() << "AsagiLite:" << path_ << "queried outside its volume for" << clamped << "of"
               << count << "points.";
  }
}

std::optional<std::pair<double, double>> AsagiLiteGrid::timeWindow() const {
  if (!timeAxis_.has_value()) {
    return std::nullopt;
  }
  const auto axis = static_cast<unsigned>(*timeAxis_);
  if (axis >= view_.dimensions) {
    return std::nullopt;
  }
  if (window_.count == 0) {
    return std::nullopt;
  }
  const double first = timeOrigin_ + static_cast<double>(window_.start) * delta_[0];
  const double last = first + static_cast<double>(window_.count - 1) * delta_[0];
  return std::make_pair(first, last);
}

// ============================ Pointer variants ============================

// coordsToIndices() and coordToIndex<> are gone. Every one of their branches
// cast a std::floor / std::ceil / std::lround result straight to std::size_t
// with no clamp, so a coordinate below the grid origin -- or a non-finite one --
// indexed far out of bounds, guarded only by an assert that is compiled out in
// release. Sampling now goes through Interpolation.cpp, which clamps first.
// Nothing in the tree called them.

// getValue() is gone with the coordinate API above: its only caller would have
// been the index path, and Interpolation.cpp reads the buffer directly through
// the ArrayView. Nothing in the tree called it.

// ============================ Allocation ============================

AsagiLiteGrid::Error AsagiLiteGrid::allocateLocal(std::size_t bytes) {
  void* p = std::malloc(bytes);
  if (p == nullptr) {
    return Error::AllocationFailed;
  }
  data_ = p;
  owner_ = Owner::Malloc;
  return Error::Success;
}

AsagiLiteGrid::Error AsagiLiteGrid::allocateShared(std::size_t bytes) {
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

} // namespace seissol::reader::datafield
