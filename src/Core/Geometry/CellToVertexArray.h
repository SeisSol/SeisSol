// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_CORE_GEOMETRY_CELLTOVERTEXARRAY_H_
#define SEISSOL_SRC_CORE_GEOMETRY_CELLTOVERTEXARRAY_H_

#include "Core/Common/Constants.h"
#include "Core/Geometry/MeshReader.h"
#include "Core/Geometry/PUMLReader.h"

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <functional>
#include <vector>

namespace seissol::initializer {

// temporary struct until we have something like a lazy vector/iterator "map" (as in on-demand,
// element-wise function application)
struct CellToVertexArray {
  using CellToVertexFunction = std::function<std::array<Eigen::Vector3d, 4>(size_t)>;
  using CellToGroupFunction = std::function<int(size_t)>;

  CellToVertexArray(size_t size,
                    CellToVertexFunction elementCoordinates,
                    CellToGroupFunction elementGroups);

  size_t size;
  CellToVertexFunction elementCoordinates;
  CellToGroupFunction elementGroups;

  static CellToVertexArray fromMeshReader(const seissol::geometry::MeshReader& meshReader);
#ifdef USE_HDF
  /// The cells of a PUML mesh before PUMLReader::getMesh(), which PUMLReader::read read with
  /// `arrays`, with their vertices in the given (canonical) order rather than in the order of the
  /// mesh file. The array refers to vertexOrders, which therefore has to outlive it.
  static CellToVertexArray
      fromPUML(const seissol::geometry::PumlMesh& mesh,
               const seissol::geometry::PumlCellArrays& arrays,
               const std::vector<seissol::geometry::VertexOrder>& vertexOrders);
  static CellToVertexArray
      fromPUML(const seissol::geometry::PumlMesh& mesh,
               const seissol::geometry::PumlCellArrays& arrays,
               std::vector<seissol::geometry::VertexOrder>&& vertexOrders) = delete;
#endif
  static CellToVertexArray
      fromVectors(const std::vector<std::array<std::array<double, 3>, 4>>& vertices,
                  const std::vector<int>& groups);
  static CellToVertexArray join(std::vector<CellToVertexArray> arrays);
  /// The cells of `array` at `indices`, in that order.
  static CellToVertexArray subset(const CellToVertexArray& array, std::vector<std::size_t> indices);
};

/// Writes the coordinates of the points [first, first + count) to xyz[i * 3 + d], with the points
/// `referencePoints` (in reference coordinates) in every cell of `cells`: point p is point
/// p % referencePoints.size() of the cell p / referencePoints.size(). The transformation of a cell
/// is set up once per range and cell rather than per point.
void cellPointCoordinates(const CellToVertexArray& cells,
                          const std::vector<std::array<double, Cell::Dim>>& referencePoints,
                          std::size_t first,
                          std::size_t count,
                          double* xyz);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_CORE_GEOMETRY_CELLTOVERTEXARRAY_H_
