// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "OutputAux.h"

#include "Common/Constants.h"
#include "Common/Iterator.h"
#include "DynamicRupture/Output/DataTypes.h"
#include "DynamicRupture/Output/Geometry.h"
#include "GeneratedCode/init.h"
#include "Geometry.h"
#include "Geometry/CellTransform.h"
#include "Geometry/FaceTransform.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshTools.h"
#include "Kernels/Precision.h"
#include "Numerical/BasisFunction.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

namespace {
double distance(const double v1[2], const double v2[2]) {
  const Eigen::Vector2d vector1(v1[0], v1[1]);
  const Eigen::Vector2d vector2(v2[0], v2[1]);
  return (vector1 - vector2).norm();
}
} // namespace

namespace seissol::dr {

int getElementVertexId(int localSideId, int localFaceVertexId) {
  return MeshTools::FACE2NODES[localSideId][localFaceVertexId];
}

ExtTriangle toExtTriangle(const geometry::FaceTransform& face) {
  const auto corners = face.vertices();
  ExtTriangle triangle{};
  for (std::size_t vertex = 0; vertex < Face::NumVertices; ++vertex) {
    for (std::size_t d = 0; d < Cell::Dim; ++d) {
      triangle.point(vertex)[d] = corners[vertex](d);
    }
  }
  return triangle;
}

ExtTriangle getReferenceTriangle(std::size_t sideIdx) {
  const geometry::ReferenceFaceMap map(sideIdx);
  ExtTriangle triangle{};
  const std::array<geometry::ReferenceFaceMap::FaceVectorT, Face::NumVertices> corners{
      geometry::ReferenceFaceMap::FaceVectorT(0.0, 0.0),
      geometry::ReferenceFaceMap::FaceVectorT(1.0, 0.0),
      geometry::ReferenceFaceMap::FaceVectorT(0.0, 1.0)};
  for (std::size_t vertex = 0; vertex < Face::NumVertices; ++vertex) {
    const auto point = map.faceToCell(corners[vertex]);
    for (std::size_t d = 0; d < Cell::Dim; ++d) {
      triangle.point(vertex)[d] = point(d);
    }
  }
  return triangle;
}

CoordinateT getMidPointTriangle(const ExtTriangle& triangle) {
  CoordinateT avgPoint{};
  const auto p0 = triangle.point(0);
  const auto p1 = triangle.point(1);
  const auto p2 = triangle.point(2);
  for (std::size_t axis = 0; axis < Cell::Dim; ++axis) {
    avgPoint[axis] = (p0[axis] + p1[axis] + p2[axis]) / 3.0;
  }
  return avgPoint;
}

CoordinateT getTrianglePointByCoords(const ExtTriangle& triangle,
                                     const std::array<double, 2>& point) {
  CoordinateT avgPoint{};
  const auto& p0 = triangle.point(0);
  const auto& p1 = triangle.point(1);
  const auto& p2 = triangle.point(2);

  // barycentric coordinates
  const auto w0 = 1 - point[0] - point[1];
  const auto w1 = point[0];
  const auto w2 = point[1];
  for (std::size_t axis = 0; axis < Cell::Dim; ++axis) {
    avgPoint[axis] = w0 * p0[axis] + w1 * p1[axis] + w2 * p2[axis];
  }
  return avgPoint;
}

CoordinateT getMidPoint(const CoordinateT& p1, const CoordinateT& p2) {
  CoordinateT midPoint{};
  for (std::size_t axis = 0; axis < Cell::Dim; ++axis) {
    midPoint[axis] = 0.5 * (p1[axis] + p2[axis]);
  }
  return midPoint;
}

TriangleQuadratureData generateTriangleQuadrature() {
  TriangleQuadratureData data{};

  // Generate triangle quadrature points and weights (Factory Method)
  const auto pointsView = init::quadpoints::view::create(init::quadpoints::Values);
  const auto weightsView = init::quadweights::view::create(init::quadweights::Values);

  auto* reshapedPoints = unsafe_reshape<2>((data.points).data());
  for (size_t i = 0; i < seissol::dr::TriangleQuadratureData::Size; ++i) {
    reshapedPoints[i][0] = seissol::multisim::multisimTranspose(pointsView, i, 0);
    reshapedPoints[i][1] = seissol::multisim::multisimTranspose(pointsView, i, 1);
    data.weights[i] = weightsView(i);
  }

  return data;
}

std::pair<int, double> getNearestFacePoint(const double targetPoint[2],
                                           const double (*facePoints)[2],
                                           std::size_t numFacePoints) {

  int nearestPoint{-1};
  double shortestDistance = std::numeric_limits<double>::max();

  for (std::size_t index = 0; index < numFacePoints; ++index) {
    const double nextPoint[2] = {facePoints[index][0], facePoints[index][1]};

    const auto currentDistance = distance(targetPoint, nextPoint);
    if (shortestDistance > currentDistance) {
      shortestDistance = currentDistance;
      nearestPoint = static_cast<int>(index);
    }
  }
  return std::make_pair(nearestPoint, shortestDistance);
}

void assignNearestGaussianPoints(Receivers& geoPoints) {
  auto quadratureData = generateTriangleQuadrature();
  const double (*trianglePoints2D)[2] = unsafe_reshape<2>(quadratureData.points.data());

  for (auto& geoPoint : geoPoints) {

    const auto targetPoint2D =
        geometry::ReferenceFaceMap(geoPoint.localFaceSideId.value())
            .cellToFace(geometry::CellTransform::VectorEigenT(geoPoint.reference.data()));

    int nearestPoint{-1};
    double shortestDistance = std::numeric_limits<double>::max();
    std::tie(nearestPoint, shortestDistance) = getNearestFacePoint(
        targetPoint2D.data(), trianglePoints2D, seissol::dr::TriangleQuadratureData::Size);
    geoPoint.nearestGpIndex = nearestPoint;
  }
}

int getClosestInternalStroudGp(int nearestGpIndex, int nPoly) {
  int i1 = ((nearestGpIndex - 1) / (nPoly + 2)) + 1;
  int j1 = (nearestGpIndex - 1) % (nPoly + 2) + 1;
  if (i1 == 1) {
    i1 = i1 + 1;
  } else if (i1 == (nPoly + 2)) {
    i1 = i1 - 1;
  }

  if (j1 == 1) {
    j1 = j1 + 1;
  } else if (j1 == (nPoly + 2)) {
    j1 = j1 - 1;
  }
  return (i1 - 1) * (nPoly + 2) + j1;
}

void projectPointToFace(CoordinateT& point,
                        const ExtTriangle& face,
                        const CoordinateT& faceNormal) {
  const auto distance = getDistanceFromPointToFace(point, face, faceNormal);
  const double faceNormalLength = MeshTools::norm(faceNormal);
  const auto adjustedDistance = distance / faceNormalLength;

  for (int i = 0; i < 3; ++i) {
    point[i] += adjustedDistance * faceNormal[i];
  }
}

double getDistanceFromPointToFace(const CoordinateT& point,
                                  const ExtTriangle& face,
                                  const CoordinateT& faceNormal) {

  CoordinateT diff{0.0, 0.0, 0.0};
  MeshTools::sub(face.point(0), point, diff);

  // Note: faceNormal may not be precisely a unit vector
  const double faceNormalLength = MeshTools::norm(faceNormal);
  return MeshTools::dot(faceNormal, diff) / faceNormalLength;
}

// (NOTE: only the sign really has a meaning; except maybe for some small tolerance)
// (reason: lack of normalization, probably)
double
    isInsideFace(const CoordinateT& point, const ExtTriangle& face, const CoordinateT& faceNormal) {

  // view the triangle as an intersection of hyperplanes

  double sidemin = std::numeric_limits<double>::max();
  for (auto [i1, i2] : seissol::common::zip(std::vector{0, 1, 2}, std::vector{1, 2, 0})) {
    const auto& p1 = face.point(i1);
    const auto& p2 = face.point(i2);
    CoordinateT sidevec{0.0, 0.0, 0.0};
    CoordinateT hypersupport{0.0, 0.0, 0.0};
    MeshTools::sub(p2, p1, sidevec);
    MeshTools::cross(faceNormal, sidevec, hypersupport);
    const auto sidevalue = MeshTools::dot(hypersupport, p1);
    const auto pointvalue = MeshTools::dot(hypersupport, point);
    const auto containvalue = pointvalue - sidevalue;
    sidemin = std::min(sidemin, containvalue);
  }
  return sidemin;
}

PlusMinusBasisFunctions getPlusMinusBasisFunctions(const CoordinateT& pointCoords,
                                                   const geometry::CellTransform& plusTransform,
                                                   const geometry::CellTransform& minusTransform) {

  Eigen::Vector3d point(pointCoords[0], pointCoords[1], pointCoords[2]);
  return getPlusMinusBasisFunctions(plusTransform.spaceToRef(point),
                                    minusTransform.spaceToRef(point));
}

PlusMinusBasisFunctions
    getPlusMinusBasisFunctions(const geometry::CellTransform::VectorEigenT& plusReference,
                               const geometry::CellTransform::VectorEigenT& minusReference) {
  auto getBasisFunctions = [](const geometry::CellTransform::VectorEigenT& referenceCoords) {
    const basisFunction::SampledBasisFunctions<real> sampler(
        ConvergenceOrder, referenceCoords[0], referenceCoords[1], referenceCoords[2]);
    return sampler.data();
  };

  PlusMinusBasisFunctions basisFunctions{};
  basisFunctions.plusSide = getBasisFunctions(plusReference);
  basisFunctions.minusSide = getBasisFunctions(minusReference);

  return basisFunctions;
}

geometry::FaceTransform::FaceVectorT
    closestPointOnFace(const geometry::FaceTransform& face,
                       const geometry::FaceTransform::VectorT& point,
                       const geometry::FaceTransform::FaceVectorT& start) {
  // converges quadratically where the point is on the face, and linearly, at a rate of its
  // distance over the radius of curvature, where it is off it
  constexpr int MaxIterations = 50;
  constexpr double Tolerance = 1e-14;
  geometry::FaceTransform::FaceVectorT chi = start;
  for (int iteration = 0; iteration < MaxIterations; ++iteration) {
    const auto jacobian = face.refToSpaceJacobian(chi);
    const geometry::FaceTransform::VectorT residual = point - face.refToSpace(chi);
    const geometry::FaceTransform::FaceVectorT step =
        (jacobian.transpose() * jacobian).ldlt().solve(jacobian.transpose() * residual);
    chi += step;
    if (step.norm() < Tolerance) {
      break;
    }
  }
  return chi;
}

real computeTriangleArea(ExtTriangle& triangle) {
  const auto p0 = Eigen::Vector3d(triangle.point(0).data());
  const auto p1 = Eigen::Vector3d(triangle.point(1).data());
  const auto p2 = Eigen::Vector3d(triangle.point(2).data());

  const auto vector1 = p1 - p0;
  const auto vector2 = p2 - p0;
  const auto normal = vector1.cross(vector2);
  return 0.5 * normal.norm();
}

std::size_t
    firstReceiverOfCell(std::size_t cell, std::size_t pointsPerCell, std::size_t simulationCount) {
  return cell * pointsPerCell * simulationCount;
}

int faultTagOfCell(const Receivers& receivers,
                   std::size_t cell,
                   std::size_t pointsPerCell,
                   std::size_t simulationCount) {
  return receivers[firstReceiverOfCell(cell, pointsPerCell, simulationCount)].faultTag;
}

std::size_t globalFaceIdOfCell(const Receivers& receivers,
                               std::size_t cell,
                               std::size_t pointsPerCell,
                               std::size_t simulationCount) {
  return receivers[firstReceiverOfCell(cell, pointsPerCell, simulationCount)].globalFaultFaceId();
}
} // namespace seissol::dr
