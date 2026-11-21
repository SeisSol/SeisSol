// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "CellToVertexArray.h"

#include "Core/Common/Constants.h"
#include "Core/Geometry/CellTransform.h"
#include "Core/Geometry/MeshDefinition.h"
#include "Core/Geometry/MeshReader.h"
#include "Core/Geometry/PUMLReader.h"

#include <Eigen/Dense>
#include <array>
#include <cassert>
#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef USE_HDF
// PUML.h needs to be included before Downward.h

#include <PUML/Downward.h>
#include <PUML/Types.h>
#endif

namespace seissol::initializer {

CellToVertexArray::CellToVertexArray(size_t size,
                                     CellToVertexFunction elementCoordinates,
                                     CellToGroupFunction elementGroups)
    : size(size), elementCoordinates(std::move(elementCoordinates)),
      elementGroups(std::move(elementGroups)) {}

CellToVertexArray
    CellToVertexArray::fromMeshReader(const seissol::geometry::MeshReader& meshReader) {
  const auto& elements = meshReader.getElements();
  const auto& vertices = meshReader.getVertices();

  return CellToVertexArray(
      elements.size(),
      [&](size_t index) {
        std::array<Eigen::Vector3d, 4> verts;
        for (size_t i = 0; i < Cell::NumVertices; ++i) {
          const auto vindex = elements[index].vertices[i];
          const auto& vertex = vertices[vindex];
          verts[i] << vertex.coords[0], vertex.coords[1], vertex.coords[2];
        }
        return verts;
      },
      [&](size_t index) { return elements[index].group; });
}

#ifdef USE_HDF
CellToVertexArray
    CellToVertexArray::fromPUML(const seissol::geometry::PumlMesh& mesh,
                                const seissol::geometry::PumlCellArrays& arrays,
                                const std::vector<seissol::geometry::VertexOrder>& vertexOrders) {
  const auto groups = mesh.data(arrays.group);
  const auto& elements = mesh.cells();
  const auto& vertices = mesh.vertices();
  assert(vertexOrders.size() == elements.size());
  return CellToVertexArray(
      elements.size(),
      [&](size_t cell) {
        std::array<Eigen::Vector3d, 4> x;
        PUML::LocalId vertLids[Cell::NumVertices]{};
        PUML::Downward::vertices(mesh, elements[cell], vertLids);
        const auto& order = vertexOrders[cell];
        for (std::size_t vtx = 0; vtx < Cell::NumVertices; ++vtx) {
          for (std::size_t d = 0; d < Cell::Dim; ++d) {
            x[vtx](d) = vertices[vertLids[order[vtx]]].coordinate()[d];
          }
        }
        return x;
      },
      [groups](size_t cell) { return groups[cell]; });
}
#endif

CellToVertexArray CellToVertexArray::fromVectors(
    const std::vector<std::array<std::array<double, Cell::Dim>, Cell::NumVertices>>& vertices,
    const std::vector<int>& groups) {
  assert(vertices.size() == groups.size());

  return CellToVertexArray(
      vertices.size(),
      [&](size_t idx) {
        std::array<Eigen::Vector3d, Cell::NumVertices> verts;
        for (size_t i = 0; i < Cell::NumVertices; ++i) {
          verts[i] << vertices[idx][i][0], vertices[idx][i][1], vertices[idx][i][2];
        }
        return verts;
      },
      [&](size_t i) { return groups[i]; });
}

CellToVertexArray CellToVertexArray::join(std::vector<CellToVertexArray> arrays) {
  std::size_t totalSize = 0;
  std::vector<std::size_t> sizes(arrays.size());
  std::vector<std::size_t> offsets(arrays.size());
  for (std::size_t i = 0; i < arrays.size(); ++i) {
    offsets[i] = totalSize;
    totalSize += arrays[i].size;
    sizes[i] = totalSize;
  }
  return CellToVertexArray(
      totalSize,
      [=](size_t idx) {
        for (std::size_t i = 0; i < sizes.size(); ++i) {
          if (idx < sizes[i]) {
            return arrays[i].elementCoordinates(idx - offsets[i]);
          }
        }
        throw std::out_of_range(std::to_string(idx) + " vs " + std::to_string(totalSize));
      },
      [=](size_t idx) {
        for (std::size_t i = 0; i < sizes.size(); ++i) {
          if (idx < sizes[i]) {
            return arrays[i].elementGroups(idx - offsets[i]);
          }
        }
        throw std::out_of_range(std::to_string(idx) + " vs " + std::to_string(totalSize));
      });
}

CellToVertexArray CellToVertexArray::subset(const CellToVertexArray& array,
                                            std::vector<std::size_t> indices) {
  const auto shared = std::make_shared<const std::vector<std::size_t>>(std::move(indices));
  return CellToVertexArray(
      shared->size(),
      [array, shared](size_t idx) { return array.elementCoordinates((*shared)[idx]); },
      [array, shared](size_t idx) { return array.elementGroups((*shared)[idx]); });
}

void cellPointCoordinates(const CellToVertexArray& cells,
                          const std::vector<std::array<double, Cell::Dim>>& referencePoints,
                          std::size_t first,
                          std::size_t count,
                          double* xyz) {
  const auto pointsPerCell = referencePoints.size();
  std::optional<seissol::geometry::AffineTransform> transform;
  std::size_t transformCell = 0;
  for (std::size_t i = 0; i < count; ++i) {
    const auto point = first + i;
    const auto cell = point / pointsPerCell;
    if (!transform.has_value() || transformCell != cell) {
      transform.emplace(cells.elementCoordinates(cell));
      transformCell = cell;
    }
    const auto transformed = transform->refToSpace(referencePoints[point % pointsPerCell]);
    for (std::size_t d = 0; d < Cell::Dim; ++d) {
      xyz[i * Cell::Dim + d] = transformed[d];
    }
  }
}

} // namespace seissol::initializer
