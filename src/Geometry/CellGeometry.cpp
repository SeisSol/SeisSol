// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "CellGeometry.h"

#include "Geometry/CellTransform.h"
#include "Geometry/FaceTransform.h"
#include "Geometry/IsoparametricTransform.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <utils/logger.h>
#include <vector>

namespace seissol::geometry {

namespace {

auto isoparametricOf(std::size_t id, const MeshReader& mesh) -> IsoparametricTransform {
  const auto nodes = mesh.cellNodes(id);
  std::vector<CellTransform::VectorEigenT> eigenNodes;
  eigenNodes.reserve(nodes.size());
  for (const auto& node : nodes) {
    eigenNodes.emplace_back(node[0], node[1], node[2]);
  }
  return {mesh.geometryOrder(), eigenNodes};
}

} // namespace

auto cellTransformOf(std::size_t id, const MeshReader& mesh) -> std::unique_ptr<CellTransform> {
  if (mesh.geometryOrder() <= 1) {
    return std::make_unique<AffineTransform>(AffineTransform::fromMeshCell(id, mesh));
  }
  return std::make_unique<IsoparametricTransform>(isoparametricOf(id, mesh));
}

auto faceTransformOf(std::size_t id,
                     std::size_t side,
                     const MeshReader& mesh,
                     FaceOrientation orientation) -> std::unique_ptr<FaceTransform> {
  if (mesh.geometryOrder() <= 1) {
    return std::make_unique<AffineFaceTransform>(
        AffineFaceTransform::fromMeshCell(id, side, mesh, orientation));
  }
  return std::make_unique<IsoparametricFaceTransform>(isoparametricOf(id, mesh),
                                                      ReferenceFaceMap(side, orientation));
}

auto faultFaceTransformOf(std::size_t faultId, const MeshReader& mesh)
    -> std::unique_ptr<FaceTransform> {
  const auto& fault = mesh.getFault().at(faultId);
  if (fault.element.hasValue()) {
    return faceTransformOf(fault.element.value(), fault.side, mesh);
  }
  if (!fault.neighborElement.hasValue()) {
    logError() << "Fault face" << faultId << "has no adjacent cell on this rank.";
  }
  const auto neighbor = fault.neighborElement.value();
  const auto orientation = static_cast<FaceOrientation>(
      mesh.getElements()[neighbor].sideOrientations[fault.neighborSide]);
  return faceTransformOf(neighbor, fault.neighborSide, mesh, orientation);
}

auto faultFrameAt(const FaceTransform& face,
                  const Fault& fault,
                  const FaceTransform::FaceVectorT& point) -> FaultFrame {
  const auto toEigen = [](const CoordinateT& vector) {
    return CellTransform::VectorEigenT(vector[0], vector[1], vector[2]);
  };
  const CellTransform::VectorEigenT normal = face.normal(point);
  const double surfaceJacobian = normal.norm();
  if (dynamic_cast<const AffineFaceTransform*>(&face) != nullptr) {
    // to the bit what a straight-sided mesh has
    return {
        toEigen(fault.normal), toEigen(fault.tangent1), toEigen(fault.tangent2), surfaceJacobian};
  }
  const CellTransform::VectorEigenT unitNormal = normal / surfaceJacobian;
  if (unitNormal.dot(toEigen(fault.normal)) <= 0) {
    logError() << "A fault face turns over: its normal at a point points against the one of the "
                  "plane through its vertices.";
  }
  const auto straightTangent = toEigen(fault.tangent1);
  const CellTransform::VectorEigenT tangent1 =
      (straightTangent - straightTangent.dot(unitNormal) * unitNormal).normalized();
  return {unitNormal, tangent1, unitNormal.cross(tangent1), surfaceJacobian};
}

auto relativeThickness(const CellTransform& transform,
                       const std::array<CellTransform::VectorEigenT, Cell::NumVertices>& vertices)
    -> double {
  if (dynamic_cast<const AffineTransform*>(&transform) != nullptr) {
    return 1.0;
  }
  const auto smallestSingularValue = [](const CellTransform::MatrixEigenT& jacobian) {
    return Eigen::JacobiSVD<CellTransform::MatrixEigenT>(jacobian).singularValues().minCoeff();
  };
  const double straight = smallestSingularValue(AffineTransform(vertices).jacobian());
  // the lattice of order four holds the nodes of a cell of order two and the points between them
  auto points = IsoparametricTransform::latticeNodes(4);
  points.emplace_back(Cell::ReferenceBarycenter.data());
  double thickness = 1.0;
  for (const auto& point : points) {
    thickness =
        std::min(thickness, smallestSingularValue(transform.refToSpaceJacobian(point)) / straight);
  }
  return thickness;
}

} // namespace seissol::geometry
