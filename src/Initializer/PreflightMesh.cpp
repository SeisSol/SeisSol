// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "PreflightMesh.h"

#include "Common/Constants.h"
#include "DynamicRupture/Output/Geometry.h"
#include "DynamicRupture/Output/OutputAux.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshTools.h"
#include "Geometry/PUMLReader.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/InitProcedure/Internal/FaceTypeCheck.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "Initializer/PointMapper.h"
#include "Parallel/OpenMP.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <hdf5.h>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer::meshscan {

namespace {

using parameters::BoundaryFormat;

//! @brief A dataset of the mesh file, which the scan reads a run of rows of at a time.
class Dataset {
  public:
  Dataset(hid_t file, const char* name) : name_(name) {
    if (H5Lexists(file, name, H5P_DEFAULT) <= 0) {
      return;
    }
    dataset_ = H5Dopen2(file, name, H5P_DEFAULT);
    if (dataset_ < 0) {
      return;
    }
    const hid_t space = H5Dget_space(dataset_);
    const int rank = H5Sget_simple_extent_ndims(space);
    if (rank > 0) {
      dims_.resize(static_cast<std::size_t>(rank));
      H5Sget_simple_extent_dims(space, dims_.data(), nullptr);
    }
    H5Sclose(space);
  }
  ~Dataset() {
    if (dataset_ >= 0) {
      H5Dclose(dataset_);
    }
  }
  Dataset(const Dataset&) = delete;
  Dataset& operator=(const Dataset&) = delete;
  Dataset(Dataset&&) = delete;
  Dataset& operator=(Dataset&&) = delete;

  [[nodiscard]] bool valid() const { return dataset_ >= 0 && !dims_.empty(); }
  [[nodiscard]] std::size_t rows() const { return dims_.empty() ? 0 : dims_[0]; }
  [[nodiscard]] std::size_t rank() const { return dims_.size(); }
  [[nodiscard]] std::size_t columns() const {
    std::size_t columns = 1;
    for (std::size_t i = 1; i < dims_.size(); ++i) {
      columns *= dims_[i];
    }
    return columns;
  }

  //! @brief Whether the dataset holds signed integers, as PUMLReader::read asks it.
  [[nodiscard]] bool signedIntegers() const {
    const hid_t type = H5Dget_type(dataset_);
    const bool isSigned = H5Tget_class(type) == H5T_INTEGER && H5Tget_sign(type) == H5T_SGN_2;
    H5Tclose(type);
    return isSigned;
  }

  //! @brief Reads the rows [first, first + count) into @p data , converted to @p memoryType .
  void read(std::size_t first, std::size_t count, hid_t memoryType, void* data) const {
    if (count == 0) {
      return;
    }
    const hid_t fileSpace = H5Dget_space(dataset_);
    std::vector<hsize_t> start(dims_.size(), 0);
    std::vector<hsize_t> extent = dims_;
    start[0] = first;
    extent[0] = count;
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, start.data(), nullptr, extent.data(), nullptr);
    const hid_t memorySpace =
        H5Screate_simple(static_cast<int>(extent.size()), extent.data(), nullptr);
    const auto status = H5Dread(dataset_, memoryType, memorySpace, fileSpace, H5P_DEFAULT, data);
    H5Sclose(memorySpace);
    H5Sclose(fileSpace);
    if (status < 0) {
      logError() << "Could not read the rows" << first << "to" << first + count << "of" << name_;
    }
  }

  private:
  std::string name_;
  hid_t dataset_{-1};
  std::vector<hsize_t> dims_;
};

/**
 * @brief The points in a regular grid over their bounding box, of about one point per cell along
 * the axes the points spread over, to find those in a box quickly.
 */
class PointGrid {
  public:
  explicit PointGrid(const std::vector<Eigen::Vector3d>& points) {
    if (points.empty()) {
      return;
    }
    lo_ = points.front();
    hi_ = points.front();
    for (const auto& point : points) {
      lo_ = lo_.cwiseMin(point);
      hi_ = hi_.cwiseMax(point);
    }
    const Eigen::Vector3d extent = hi_ - lo_;
    double volume = 1;
    int spread = 0;
    for (int i = 0; i < 3; ++i) {
      if (extent(i) > 0) {
        volume *= extent(i);
        ++spread;
      }
    }
    const double side =
        spread > 0 ? std::pow(volume / static_cast<double>(points.size()), 1.0 / spread) : 1.0;
    constexpr std::size_t MaxCells = 1024;
    for (int i = 0; i < 3; ++i) {
      dims_[i] = extent(i) > 0
                     ? std::clamp<std::size_t>(
                           static_cast<std::size_t>(std::ceil(extent(i) / side)), 1, MaxCells)
                     : 1;
      cell_(i) = extent(i) > 0 ? extent(i) / static_cast<double>(dims_[i]) : 1.0;
    }

    std::vector<std::size_t> cellOf(points.size());
    start_.assign(dims_[0] * dims_[1] * dims_[2] + 1, 0);
    for (std::size_t i = 0; i < points.size(); ++i) {
      cellOf[i] = flat(index(points[i](0), 0), index(points[i](1), 1), index(points[i](2), 2));
      ++start_[cellOf[i] + 1];
    }
    std::partial_sum(start_.begin(), start_.end(), start_.begin());
    order_.resize(points.size());
    auto next = start_;
    for (std::size_t i = 0; i < points.size(); ++i) {
      order_[next[cellOf[i]]++] = i;
    }
  }

  //! @brief Calls @p visit with the index of every point in a grid cell that the box [a, b] meets.
  template <typename F>
  void forEachNear(const Eigen::Vector3d& a, const Eigen::Vector3d& b, F&& visit) const {
    if (order_.empty()) {
      return;
    }
    for (int i = 0; i < 3; ++i) {
      if (b(i) < lo_(i) || a(i) > hi_(i)) {
        return;
      }
    }
    const std::array<std::size_t, 3> from{index(a(0), 0), index(a(1), 1), index(a(2), 2)};
    const std::array<std::size_t, 3> to{index(b(0), 0), index(b(1), 1), index(b(2), 2)};
    for (std::size_t x = from[0]; x <= to[0]; ++x) {
      for (std::size_t y = from[1]; y <= to[1]; ++y) {
        for (std::size_t z = from[2]; z <= to[2]; ++z) {
          const auto cell = flat(x, y, z);
          for (auto k = start_[cell]; k < start_[cell + 1]; ++k) {
            std::forward<F>(visit)(order_[k]);
          }
        }
      }
    }
  }

  private:
  [[nodiscard]] std::size_t index(double value, int axis) const {
    const double position = std::floor((value - lo_(axis)) / cell_(axis));
    if (!(position > 0)) {
      return 0;
    }
    if (!(position < static_cast<double>(dims_[axis]))) {
      return dims_[axis] - 1;
    }
    return static_cast<std::size_t>(position);
  }
  [[nodiscard]] std::size_t flat(std::size_t x, std::size_t y, std::size_t z) const {
    return (x * dims_[1] + y) * dims_[2] + z;
  }

  Eigen::Vector3d lo_{Eigen::Vector3d::Zero()};
  Eigen::Vector3d hi_{Eigen::Vector3d::Zero()};
  Eigen::Vector3d cell_{Eigen::Vector3d::Ones()};
  std::array<std::size_t, 3> dims_{1, 1, 1};
  std::vector<std::size_t> start_;
  std::vector<std::size_t> order_;
};

//! @brief Where the run would take a point to lie: the key of findUniqueMeshIds and the cell.
struct Candidate {
  double key{std::numeric_limits<double>::infinity()};
  std::size_t cell{std::numeric_limits<std::size_t>::max()};

  [[nodiscard]] bool betterThan(const Candidate& other) const {
    return key < other.key || (key == other.key && cell < other.cell);
  }
};

//! @brief A face of dynamic rupture, as one of its cells sees it.
struct FaultFace {
  //! the vertices of the face, in the order of the face (MeshTools::FACE2NODES) in its cell
  std::array<std::uint64_t, 3> vertices{};
  std::size_t cell{0};
  bool plus{false};

  //! the vertices, sorted, which both cells of the face agree on
  [[nodiscard]] std::array<std::uint64_t, 3> key() const {
    auto sorted = vertices;
    std::ranges::sort(sorted);
    return sorted;
  }
};

//! @brief What a thread gathers over the cells it scans.
struct Gathered {
  std::array<double, 3> min{std::numeric_limits<double>::infinity(),
                            std::numeric_limits<double>::infinity(),
                            std::numeric_limits<double>::infinity()};
  std::array<double, 3> max{-std::numeric_limits<double>::infinity(),
                            -std::numeric_limits<double>::infinity(),
                            -std::numeric_limits<double>::infinity()};
  std::map<std::uint32_t, std::size_t> tags;
  std::map<std::uint32_t, std::size_t> unknownTags;
  std::map<FaceType, std::size_t> faceTypes;
  std::map<int, std::size_t> groups;
  std::map<int, std::uint32_t> groupFaceTypes;
  std::vector<FaultFace> faultFaces;
  std::size_t invalidCells{0};
  std::size_t misorientedCells{0};
  std::optional<std::size_t> firstMisorientedCell;
  double smallestInsphere{std::numeric_limits<double>::infinity()};
  std::size_t smallestCell{0};
  double largestInsphere{0};
};

/**
 * @brief The side a fault face has in the run, from the cell with the vertices @p cell , in the
 * canonical order, of which it is the face @p face (MeshTools numbering): as
 * MeshReader::extractFaultInformation tells it.
 */
bool isPlusSide(const std::array<CoordinateT, Cell::NumVertices>& cell,
                std::size_t face,
                const MeshScanSettings& settings) {
  const auto& nodes = MeshTools::FACE2NODES[face];
  CoordinateT chiVec{};
  CoordinateT tauVec{};
  MeshTools::sub(cell[nodes[1]], cell[nodes[0]], chiVec);
  MeshTools::sub(cell[nodes[2]], cell[nodes[0]], tauVec);
  CoordinateT normal{};
  MeshTools::cross(chiVec, tauVec, normal);
  MeshTools::mul(normal, 1.0 / MeshTools::norm(normal), normal);

  const CoordinateT refPoint{
      settings.referencePoint(0), settings.referencePoint(1), settings.referencePoint(2)};
  if (settings.refPointMethod == parameters::RefPointMethod::Point) {
    CoordinateT tmp1{};
    CoordinateT tmp2{};
    MeshTools::sub(refPoint, cell[nodes[0]], tmp1);
    MeshTools::sub(cell[MeshTools::FACE2MISSINGNODE[face]], cell[nodes[0]], tmp2);
    return MeshTools::dot(tmp1, normal) * MeshTools::dot(tmp2, normal) > 0;
  }
  return MeshTools::dot(refPoint, normal) > 0;
}

//! @brief The point of the file vertex @p vertex , displaced and scaled as the run does.
CoordinateT transformed(const double* vertex, const MeshScanSettings& settings) {
  // as MeshReader::displaceMesh and MeshReader::scaleMesh, one after the other
  CoordinateT displaced{};
  for (std::size_t i = 0; i < Cell::Dim; ++i) {
    displaced[i] = vertex[i] + settings.displacement[i];
  }
  Eigen::Vector3d point;
  point << displaced[0], displaced[1], displaced[2];
  const Eigen::Vector3d result = settings.scaling * point;
  return {result[0], result[1], result[2]};
}

/**
 * @brief The height of the triangle @p a , @p b , @p c at the horizontal position of @p point , if
 * the vertical through the point meets the triangle.
 */
std::optional<double> heightAt(const Eigen::Vector3d& point,
                               const CoordinateT& a,
                               const CoordinateT& b,
                               const CoordinateT& c) {
  const double det = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
  const double scale = std::max(
      {std::abs(b[0] - a[0]), std::abs(c[0] - a[0]), std::abs(b[1] - a[1]), std::abs(c[1] - a[1])});
  // a face that stands upright has no height to give
  if (!(std::abs(det) > 1e-12 * scale * scale)) {
    return {};
  }
  const double s = ((point(0) - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (point(1) - a[1])) / det;
  const double t = ((b[0] - a[0]) * (point(1) - a[1]) - (point(0) - a[0]) * (b[1] - a[1])) / det;
  constexpr double Slack = 1e-9;
  if (s < -Slack || t < -Slack || s + t > 1 + Slack) {
    return {};
  }
  return a[2] + s * (b[2] - a[2]) + t * (c[2] - a[2]);
}

//! @brief Where the run puts the on-fault receivers, among the deduplicated @p faces .
void locateFaultPoints(FaultPoints& faultPoints,
                       const std::vector<FaultFace>& faces,
                       const std::vector<double>& fileVertices,
                       const MeshScanSettings& settings) {
  faultPoints.found.assign(faultPoints.points.size(), std::nullopt);
  if (faces.empty() || faultPoints.points.empty()) {
    return;
  }

  // the faces as the run sees them: from the plus side, with its normal
  struct Face {
    dr::ExtTriangle triangle;
    CoordinateT normal{};
    std::array<double, 3> lo{};
    std::array<double, 3> hi{};
    std::size_t cell{0};
    double size{0};
  };
  std::vector<Face> prepared(faces.size());
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < faces.size(); ++i) {
    const auto& face = faces[i];
    std::array<CoordinateT, 3> corners{};
    for (std::size_t k = 0; k < 3; ++k) {
      corners[k] = transformed(&fileVertices[3 * face.vertices[k]], settings);
    }
    // a face only the minus side of which is known has the normal of the other side, which
    // points into the cell (see MeshReader::extractFaultInformation); its corners are turned
    // around with it
    if (!face.plus) {
      std::swap(corners[1], corners[2]);
    }
    CoordinateT ab{};
    CoordinateT ac{};
    MeshTools::sub(corners[1], corners[0], ab);
    MeshTools::sub(corners[2], corners[0], ac);
    CoordinateT normal{};
    MeshTools::cross(ab, ac, normal);
    const double length = MeshTools::norm(normal);
    MeshTools::mul(normal, 1.0 / length, normal);

    auto& result = prepared[i];
    result.triangle = dr::ExtTriangle(corners[0], corners[1], corners[2]);
    result.normal = normal;
    result.cell = face.cell;
    result.size = std::sqrt(0.5 * length);
    for (std::size_t d = 0; d < 3; ++d) {
      result.lo[d] = std::min({corners[0][d], corners[1][d], corners[2][d]});
      result.hi[d] = std::max({corners[0][d], corners[1][d], corners[2][d]});
    }
  }

  // The run takes the closest face whose prism holds the point. A face can only be closer than
  // the closest one so far where its bounding box is, which rules out most faces once a close one
  // is known; the faces near the point come first, in the order of their centers along each axis.
  std::array<std::vector<std::size_t>, 3> byAxis;
  for (std::size_t d = 0; d < 3; ++d) {
    byAxis[d].resize(prepared.size());
    std::iota(byAxis[d].begin(), byAxis[d].end(), 0);
    std::ranges::sort(byAxis[d], [&](std::size_t a, std::size_t b) {
      return prepared[a].lo[d] + prepared[a].hi[d] < prepared[b].lo[d] + prepared[b].hi[d];
    });
  }

#pragma omp parallel for schedule(dynamic)
  for (std::size_t p = 0; p < faultPoints.points.size(); ++p) {
    const auto& point = faultPoints.points[p];
    double best = std::numeric_limits<double>::infinity();
    std::optional<std::size_t> bestFace;

    const auto test = [&](std::size_t f) {
      const auto& face = prepared[f];
      double gap = 0;
      for (std::size_t d = 0; d < 3; ++d) {
        const double outside = std::max({face.lo[d] - point[d], point[d] - face.hi[d], 0.0});
        gap += outside * outside;
      }
      if (std::sqrt(gap) >= best) {
        return;
      }
      // as PickPointBuilder::findClosestFaultIndex
      if (dr::isInsideFace(point, face.triangle, face.normal) > -1e-12) {
        const double distance =
            std::abs(dr::getDistanceFromPointToFace(point, face.triangle, face.normal));
        if (distance < best || (distance == best && bestFace.has_value() && f < *bestFace)) {
          best = distance;
          bestFace = f;
        }
      }
    };

    // the faces around the point along each axis first, then all of them
    for (std::size_t d = 0; d < 3; ++d) {
      const auto& order = byAxis[d];
      const auto position = std::ranges::lower_bound(order, 2 * point[d], {}, [&](std::size_t f) {
        return prepared[f].lo[d] + prepared[f].hi[d];
      });
      const auto center = static_cast<std::size_t>(position - order.begin());
      constexpr std::size_t Window = 64;
      const auto from = center > Window ? center - Window : 0;
      const auto to = std::min(order.size(), center + Window);
      for (auto k = from; k < to; ++k) {
        test(order[k]);
      }
    }
    for (std::size_t f = 0; f < prepared.size(); ++f) {
      test(f);
    }

    if (bestFace.has_value()) {
      const auto& face = prepared[*bestFace];
      faultPoints.found[p] =
          FaultPoints::Found{.cell = face.cell, .distance = best, .faceSize = face.size};
    }
  }
}

} // namespace

MeshScan scanMesh(const MeshScanSettings& settings,
                  CellPoints& receivers,
                  CellPoints& sources,
                  FaultPoints& faultReceivers) {
  const auto started = std::chrono::steady_clock::now();
  MeshScan scan;

  // HDF5 would print its error stack for every dataset that is not there
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  const hid_t file = H5Fopen(settings.file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    logError() << "Could not open the mesh file" << settings.file;
  }

  scan.boundaryFormat = settings.boundaryFormat == BoundaryFormat::Auto
                            ? geometry::boundaryFormatOf(file)
                            : settings.boundaryFormat;

  // the vertices, as the file gives them; the cells are displaced and scaled as they are read
  std::vector<double> fileVertices;
  {
    const Dataset dataset(file, "/geometry");
    if (!dataset.valid() || dataset.columns() != Cell::Dim) {
      logError() << "The mesh file" << settings.file << "has no vertices (/geometry) to scan.";
    }
    scan.vertices = dataset.rows();
    fileVertices.resize(Cell::Dim * scan.vertices);
    dataset.read(0, scan.vertices, H5T_NATIVE_DOUBLE, fileVertices.data());
  }

  const Dataset connect(file, "/connect");
  const Dataset group(file, "/group");
  const Dataset boundary(file, "/boundary");
  if (!connect.valid() || connect.columns() != Cell::NumVertices) {
    logError() << "The mesh file" << settings.file << "has no cells (/connect) to scan.";
  }
  scan.cells = connect.rows();
  const bool hasGroups = group.valid() && group.rows() == scan.cells;
  const auto boundaryValuesPerCell =
      scan.boundaryFormat == BoundaryFormat::I32x4 ? std::size_t{Cell::NumFaces} : std::size_t{1};
  const bool hasBoundaries = boundary.valid() && boundary.rows() == scan.cells &&
                             boundary.columns() == boundaryValuesPerCell;
  // as PUMLReader::read: the values with the bits they have in the file
  const hid_t boundaryType =
      hasBoundaries && boundary.signedIntegers() ? H5T_NATIVE_INT64 : H5T_NATIVE_UINT64;

  for (auto* points : {&receivers, &sources}) {
    points->found.assign(points->points.size(), CellPoints::Found{});
  }
  const PointGrid receiverGrid(receivers.points);
  const PointGrid sourceGrid(sources.points);
  std::vector<Candidate> receiverBest(receivers.points.size());
  std::vector<Candidate> sourceBest(sources.points.size());
  std::vector<CellPlanes::Point> receivers1;
  std::vector<CellPlanes::Point> sources1;
  for (const auto& point : receivers.points) {
    receivers1.push_back(CellPlanes::homogeneous(point));
  }
  for (const auto& point : sources.points) {
    sources1.push_back(CellPlanes::homogeneous(point));
  }

  std::vector<Gathered> gathered(OpenMP::threadCount());

  const auto chunk = std::max<std::size_t>(settings.chunk, 1);
  std::vector<std::uint64_t> cells;
  std::vector<int> groups;
  std::vector<std::uint64_t> boundaries;
  CellChunk visited;
  for (std::size_t first = 0; first < scan.cells; first += chunk) {
    const auto count = std::min(chunk, scan.cells - first);
    if (settings.visitCells) {
      visited.first = first;
      visited.vertices.assign(count, {});
      visited.valid.assign(count, 0);
    }
    cells.resize(Cell::NumVertices * count);
    connect.read(first, count, H5T_NATIVE_UINT64, cells.data());
    groups.assign(count, 0);
    if (hasGroups) {
      group.read(first, count, H5T_NATIVE_INT, groups.data());
    }
    boundaries.assign(boundaryValuesPerCell * count, 0);
    if (hasBoundaries) {
      boundary.read(first, count, boundaryType, boundaries.data());
    }

#pragma omp parallel for schedule(static)
    for (std::size_t local = 0; local < count; ++local) {
      auto& mine = gathered[OpenMP::threadId()];
      const auto cell = first + local;

      std::array<std::uint64_t, Cell::NumVertices> ids{};
      bool valid = true;
      for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
        ids[k] = cells[Cell::NumVertices * local + k];
        valid = valid && ids[k] < scan.vertices;
      }
      if (!valid) {
        ++mine.invalidCells;
        continue;
      }

      // the vertices in the order the run puts them in, displaced and scaled as there
      std::array<const double*, Cell::NumVertices> fileCoords{};
      for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
        fileCoords[k] = &fileVertices[Cell::Dim * ids[k]];
      }
      const auto order = geometry::canonicalVertexOrder(ids, fileCoords);
      std::array<CoordinateT, Cell::NumVertices> vertices{};
      std::array<std::uint64_t, Cell::NumVertices> vertexIds{};
      for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
        vertices[k] = transformed(fileCoords[order[k]], settings);
        vertexIds[k] = ids[order[k]];
      }
      if (settings.visitCells) {
        visited.vertices[local] = vertices;
        visited.valid[local] = 1;
      }

      // as MeshReader::verifyMeshOrientation
      Eigen::Matrix<double, 4, 4> orientation;
      for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
        orientation.row(static_cast<Eigen::Index>(k)) << vertices[k][0], vertices[k][1],
            vertices[k][2], 1;
      }
      if (orientation.determinant() >= 0) {
        ++mine.misorientedCells;
        if (!mine.firstMisorientedCell.has_value() || cell < *mine.firstMisorientedCell) {
          mine.firstMisorientedCell = cell;
        }
      }

      std::array<double, 3> lo{};
      std::array<double, 3> hi{};
      for (std::size_t d = 0; d < 3; ++d) {
        lo[d] = std::min({vertices[0][d], vertices[1][d], vertices[2][d], vertices[3][d]});
        hi[d] = std::max({vertices[0][d], vertices[1][d], vertices[2][d], vertices[3][d]});
        mine.min[d] = std::min(mine.min[d], lo[d]);
        mine.max[d] = std::max(mine.max[d], hi[d]);
      }

      std::array<const CoordinateT*, Cell::NumVertices> corners{};
      for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
        corners[k] = &vertices[k];
      }
      const CellPlanes planes(corners);

      // the size of the cell: the radius of its insphere, three times its volume by its surface;
      // and the shortest normal of a face, which is twice as long as the face is large
      double surface = 0;
      double shortestNormal = std::numeric_limits<double>::infinity();
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        const auto& nodes = MeshTools::FACE2NODES[face];
        CoordinateT ab{};
        CoordinateT ac{};
        CoordinateT n{};
        MeshTools::sub(vertices[nodes[1]], vertices[nodes[0]], ab);
        MeshTools::sub(vertices[nodes[2]], vertices[nodes[0]], ac);
        MeshTools::cross(ab, ac, n);
        surface += 0.5 * MeshTools::norm(n);
        shortestNormal = std::min(shortestNormal, MeshTools::norm(n));
      }
      const double volume = std::abs(orientation.determinant()) / 6.0;
      const double insphere = surface > 0 ? 3.0 * volume / surface : 0.0;
      if (insphere < mine.smallestInsphere) {
        mine.smallestInsphere = insphere;
        mine.smallestCell = cell;
      }
      mine.largestInsphere = std::max(mine.largestInsphere, insphere);

      // the points in the cell, by the test of the run; the bounding box of the cell holds every
      // point that the cell takes within the tolerance, widened by how far the tolerance reaches
      const auto locate = [&](const PointGrid& grid,
                              const std::vector<CellPlanes::Point>& points1,
                              double tolerance,
                              std::vector<Candidate>& best) {
        if (points1.empty()) {
          return;
        }
        double extent = 0;
        for (std::size_t d = 0; d < 3; ++d) {
          extent = std::max(extent, hi[d] - lo[d]);
        }
        // A point that the cell takes within the tolerance lies at most tolerance / |n| beyond the
        // plane of a face with the normal n, i.e. in the cell grown by that much at each face. That
        // lies in the cell scaled about the center of its insphere by (r + reach) / r, where r is
        // the radius of the insphere and reach the farthest of these distances; its vertices move
        // out by reach / r times their distance from the center at most. Beyond that, round-off.
        const double reach = tolerance / shortestNormal;
        double margin = 1e-9 * extent + (tolerance > 0 ? reach / insphere * extent : 0.0);
        if (std::isnan(margin)) {
          // a cell without volume; the run rejects it (see MeshReader::verifyMeshOrientation)
          margin = std::numeric_limits<double>::infinity();
        }
        const Eigen::Vector3d a(lo[0] - margin, lo[1] - margin, lo[2] - margin);
        const Eigen::Vector3d b(hi[0] + margin, hi[1] + margin, hi[2] + margin);
        grid.forEachNear(a, b, [&](std::size_t p) {
          const auto fit = planes.fit(points1[p]);
          if (fit.value <= tolerance) {
            const Candidate candidate{.key = planes.key(fit), .cell = cell};
#pragma omp critical(seissol_preflight_points)
            {
              if (candidate.betterThan(best[p])) {
                best[p] = candidate;
              }
            }
          }
        });
      };
      locate(receiverGrid, receivers1, receivers.tolerance, receiverBest);
      locate(sourceGrid, sources1, sources.tolerance, sourceBest);

      // the faces
      const int cellGroup = groups[local];
      ++mine.groups[cellGroup];
      auto& groupTypes = mine.groupFaceTypes[cellGroup];
      const auto faceMap = geometry::canonicalFaceMap(order);
      for (std::size_t pumlFace = 0; pumlFace < Cell::NumFaces; ++pumlFace) {
        const auto tag = geometry::boundaryTag(
            &boundaries[boundaryValuesPerCell * local], pumlFace, scan.boundaryFormat);
        if (tag != 0) {
          ++mine.tags[tag];
        }
        const auto type = settings.faceMap.at(tag);
        if (!type.has_value()) {
          ++mine.unknownTags[tag];
          continue;
        }
        // (see MeshReader::disableDR)
        const auto faceType = *type == FaceType::DynamicRupture && !settings.dynamicRupture
                                  ? FaceType::Regular
                                  : *type;
        if (faceType != FaceType::Regular) {
          ++mine.faceTypes[faceType];
        }
        groupTypes |= internal::faceTypeBit(faceType);

        const std::size_t face = faceMap[pumlFace];
        const auto& nodes = MeshTools::FACE2NODES[face];
        if (faceType == FaceType::DynamicRupture) {
          mine.faultFaces.push_back(
              FaultFace{.vertices = {vertexIds[nodes[0]], vertexIds[nodes[1]], vertexIds[nodes[2]]},
                        .cell = cell,
                        .plus = isPlusSide(vertices, face, settings)});
        }
        if (faceType == FaceType::FreeSurface || faceType == FaceType::FreeSurfaceGravity) {
          const auto& a = vertices[nodes[0]];
          const auto& b = vertices[nodes[1]];
          const auto& c = vertices[nodes[2]];
          const Eigen::Vector3d from(std::min({a[0], b[0], c[0]}),
                                     std::min({a[1], b[1], c[1]}),
                                     -std::numeric_limits<double>::infinity());
          const Eigen::Vector3d to(std::max({a[0], b[0], c[0]}),
                                   std::max({a[1], b[1], c[1]}),
                                   std::numeric_limits<double>::infinity());
          for (auto* points : {&receivers, &sources}) {
            const auto& grid = points == &receivers ? receiverGrid : sourceGrid;
            grid.forEachNear(from, to, [&](std::size_t p) {
              const auto height = heightAt(points->points[p], a, b, c);
              if (height.has_value()) {
#pragma omp critical(seissol_preflight_surface)
                {
                  auto& surface = points->found[p].surface;
                  surface = std::max(surface.value_or(*height), *height);
                }
              }
            });
          }
        }
      }
    }
    if (settings.visitCells) {
      visited.groups.assign(groups.begin(), groups.begin() + static_cast<std::ptrdiff_t>(count));
      settings.visitCells(visited);
    }
  }
  H5Fclose(file);

  // what the threads gathered
  std::vector<FaultFace> faultFaces;
  scan.smallestInsphere = std::numeric_limits<double>::infinity();
  for (int d = 0; d < 3; ++d) {
    scan.min(d) = std::numeric_limits<double>::infinity();
    scan.max(d) = -std::numeric_limits<double>::infinity();
  }
  for (auto& part : gathered) {
    for (int d = 0; d < 3; ++d) {
      scan.min(d) = std::min(scan.min(d), part.min[d]);
      scan.max(d) = std::max(scan.max(d), part.max[d]);
    }
    for (const auto& [tag, faces] : part.tags) {
      scan.tags[tag] += faces;
    }
    for (const auto& [tag, faces] : part.unknownTags) {
      scan.unknownTags[tag] += faces;
    }
    for (const auto& [type, faces] : part.faceTypes) {
      scan.faceTypes[type] += faces;
    }
    for (const auto& [id, count] : part.groups) {
      scan.groups[id] += count;
    }
    for (const auto& [id, types] : part.groupFaceTypes) {
      scan.groupFaceTypes[id] |= types;
    }
    scan.invalidCells += part.invalidCells;
    scan.misorientedCells += part.misorientedCells;
    if (part.firstMisorientedCell.has_value() &&
        (!scan.firstMisorientedCell.has_value() ||
         *part.firstMisorientedCell < *scan.firstMisorientedCell)) {
      scan.firstMisorientedCell = part.firstMisorientedCell;
    }
    if (part.smallestInsphere < scan.smallestInsphere ||
        (part.smallestInsphere == scan.smallestInsphere && part.smallestCell < scan.smallestCell)) {
      scan.smallestInsphere = part.smallestInsphere;
      scan.smallestCell = part.smallestCell;
    }
    scan.largestInsphere = std::max(scan.largestInsphere, part.largestInsphere);
    faultFaces.insert(faultFaces.end(), part.faultFaces.begin(), part.faultFaces.end());
    part.faultFaces.clear();
    part.faultFaces.shrink_to_fit();
  }
  if (scan.cells == 0) {
    scan.min.setZero();
    scan.max.setZero();
    scan.smallestInsphere = 0;
  }

  for (std::size_t p = 0; p < receivers.points.size(); ++p) {
    if (receiverBest[p].key < std::numeric_limits<double>::infinity()) {
      receivers.found[p].cell = receiverBest[p].cell;
    }
  }
  for (std::size_t p = 0; p < sources.points.size(); ++p) {
    if (sourceBest[p].key < std::numeric_limits<double>::infinity()) {
      sources.found[p].cell = sourceBest[p].cell;
    }
  }

  // each fault face once, from its plus side where one of its cells is that (see
  // MeshReader::extractFaultInformation)
  std::ranges::sort(faultFaces, [](const FaultFace& a, const FaultFace& b) {
    const auto keyA = a.key();
    const auto keyB = b.key();
    return std::tie(keyA, b.plus, a.cell) < std::tie(keyB, a.plus, b.cell);
  });
  // a fault face has a cell on either side, of which the reference point tells one to be the plus
  // side; where it lies in the plane of the face, neither is
  for (std::size_t first = 0; first < faultFaces.size();) {
    std::size_t next = first + 1;
    while (next < faultFaces.size() && faultFaces[next].key() == faultFaces[first].key()) {
      ++next;
    }
    if (next - first == 1) {
      ++scan.oneSidedFaultFaces;
    } else if (!faultFaces[first].plus) {
      ++scan.faultFacesWithoutPlusSide;
    }
    first = next;
  }
  const auto [last, end] = std::ranges::unique(
      faultFaces, [](const FaultFace& a, const FaultFace& b) { return a.key() == b.key(); });
  faultFaces.erase(last, end);
  scan.faultFaces = faultFaces.size();
  locateFaultPoints(faultReceivers, faultFaces, fileVertices, settings);

  scan.bytes = fileVertices.size() * sizeof(double) + faultFaces.size() * sizeof(FaultFace) +
               chunk * (Cell::NumVertices * sizeof(std::uint64_t) + sizeof(int) +
                        boundaryValuesPerCell * sizeof(std::uint64_t));
  scan.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
  return scan;
}

} // namespace seissol::initializer::meshscan
