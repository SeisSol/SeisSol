#include "Geometry/MeshTools.h"
#include "Numerical_aux/BasisFunction.h"
#include "Numerical_aux/Quadrature.h"
#include "Numerical_aux/Transformation.h"
#include "OutputAux.hpp"
#include <Eigen/Dense>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <unordered_map>

namespace seissol::dr {

int getElementVertexId(int localSideId, int localFaceVertexId) {
  return MeshTools::FACE2NODES[localSideId][localFaceVertexId];
}

ExtTriangle getReferenceTriangle(int sideIdx) {
  ExtTriangle referenceFace;
  switch (sideIdx) {
  case 0:
    referenceFace.point(0) = {0.0, 0.0, 0.0};
    referenceFace.point(1) = {0.0, 1.0, 0.0};
    referenceFace.point(2) = {1.0, 0.0, 0.0};
    break;
  case 1:
    referenceFace.point(0) = {0.0, 0.0, 0.0};
    referenceFace.point(1) = {1.0, 0.0, 0.0};
    referenceFace.point(2) = {0.0, 0.0, 1.0};
    break;
  case 2:
    referenceFace.point(0) = {0.0, 0.0, 0.0};
    referenceFace.point(1) = {0.0, 0.0, 1.0};
    referenceFace.point(2) = {0.0, 1.0, 0.0};
    break;
  case 3:
    referenceFace.point(0) = {1.0, 0.0, 0.0};
    referenceFace.point(1) = {0.0, 1.0, 0.0};
    referenceFace.point(2) = {0.0, 0.0, 1.0};
    break;
  default:
    logError() << "Unknown Local Side Id. Must be 0, 1, 2 or 3";
  }

  return referenceFace;
}

ExtTriangle getGlobalTriangle(int localSideId,
                              const Element& element,
                              const std::vector<Vertex>& verticesInfo) {
  ExtTriangle triangle{};

  for (int vertexId = 0; vertexId < 3; ++vertexId) {
    const auto elementVertexId = getElementVertexId(localSideId, vertexId);
    const auto globalVertexId = element.vertices[elementVertexId];

    triangle.point(vertexId) = verticesInfo[globalVertexId].coords;
  }
  return triangle;
}

ExtVrtxCoords getMidPointTriangle(const ExtTriangle& triangle) {
  ExtVrtxCoords avgPoint{};
  const auto p0 = triangle.point(0);
  const auto p1 = triangle.point(1);
  const auto p2 = triangle.point(2);
  for (int axis = 0; axis < 3; ++axis) {
    avgPoint.coords[axis] = (p0.coords[axis] + p1.coords[axis] + p2.coords[axis]) / 3.0;
  }
  return avgPoint;
}

ExtVrtxCoords getMidPoint(const ExtVrtxCoords& p1, const ExtVrtxCoords& p2) {
  ExtVrtxCoords midPoint{};
  for (int axis = 0; axis < 3; ++axis) {
    midPoint.coords[axis] = 0.5 * (p1.coords[axis] + p2.coords[axis]);
  }
  return midPoint;
}

TriangleQuadratureData generateTriangleQuadrature(unsigned polyDegree) {
  TriangleQuadratureData data{};

  // Generate triangle quadrature points and weights (Factory Method)
  auto pointsView = init::quadpoints::view::create(const_cast<real*>(init::quadpoints::Values));
  auto weightsView = init::quadweights::view::create(const_cast<real*>(init::quadweights::Values));
  auto* reshapedPoints = unsafe_reshape<2>(&data.points[0]);
  for (size_t i = 0; i < data.size; ++i) {
    reshapedPoints[i][0] = pointsView(i, 0);
    reshapedPoints[i][1] = pointsView(i, 1);
    data.weights[i] = weightsView(i);
  }

  return data;
}

double distance(const double v1[2], const double v2[2]) {
  Eigen::Vector2d vector1(v1[0], v1[1]), vector2(v2[0], v2[1]);
  return (vector1 - vector2).norm();
}

std::pair<int, double> getNearestFacePoint(const double targetPoint[2],
                                           const double (*facePoints)[2],
                                           const unsigned numFacePoints) {

  int nearestPoint{-1};
  double shortestDistance = std::numeric_limits<double>::max();

  for (unsigned index = 0; index < numFacePoints; ++index) {
    const double nextPoint[2] = {facePoints[index][0], facePoints[index][1]};

    const auto currentDistance = distance(targetPoint, nextPoint);
    if (shortestDistance > currentDistance) {
      shortestDistance = currentDistance;
      nearestPoint = static_cast<int>(index);
    }
  }
  return std::make_pair(nearestPoint, shortestDistance);
}

void assignNearestGaussianPoints(ReceiverPoints& geoPoints) {
  auto quadratureData = generateTriangleQuadrature(CONVERGENCE_ORDER + 1);
  double(*trianglePoints2D)[2] = unsafe_reshape<2>(&quadratureData.points[0]);

  for (auto& geoPoint : geoPoints) {

    double targetPoint2D[2];
    transformations::XiEtaZeta2chiTau(
        geoPoint.localFaceSideId, geoPoint.reference.coords, targetPoint2D);

    int nearestPoint{-1};
    double shortestDistance = std::numeric_limits<double>::max();
    std::tie(nearestPoint, shortestDistance) =
        getNearestFacePoint(targetPoint2D, trianglePoints2D, quadratureData.size);
    geoPoint.nearestGpIndex = nearestPoint;
  }
}

int getClosestInternalStroudGp(int nearestGpIndex, int nPoly) {
  int i1 = int((nearestGpIndex - 1) / (nPoly + 2)) + 1;
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

void projectPointToFace(ExtVrtxCoords& point,
                        const ExtTriangle& face,
                        const VrtxCoords faceNormal) {
  const auto distance = getDistanceFromPointToFace(point, face, faceNormal);
  const double faceNormalLength = MeshTools::norm(faceNormal);
  const auto adjustedDistance = distance / faceNormalLength;

  for (int i = 0; i < 3; ++i) {
    point.coords[i] += adjustedDistance * faceNormal[i];
  }
}

double getDistanceFromPointToFace(const ExtVrtxCoords& point,
                                  const ExtTriangle& face,
                                  const VrtxCoords faceNormal) {

  VrtxCoords diff{0.0, 0.0, 0.0};
  MeshTools::sub(face.point(0).coords, point.coords, diff);

  // Note: faceNormal may not be precisely a unit vector
  const double faceNormalLength = MeshTools::norm(faceNormal);
  return MeshTools::dot(faceNormal, diff) / faceNormalLength;
}

PlusMinusBasisFunctions getPlusMinusBasisFunctions(const VrtxCoords pointCoords,
                                                   const VrtxCoords* plusElementCoords[4],
                                                   const VrtxCoords* minusElementCoords[4]) {

  Eigen::Vector3d point(pointCoords[0], pointCoords[1], pointCoords[2]);

  auto getBasisFunctions = [&point](const VrtxCoords* elementCoords[4]) {
    auto referenceCoords = transformations::tetrahedronGlobalToReference(
        *elementCoords[0], *elementCoords[1], *elementCoords[2], *elementCoords[3], point);
    basisFunction::SampledBasisFunctions<real> sampler(
        CONVERGENCE_ORDER, referenceCoords[0], referenceCoords[1], referenceCoords[2]);
    return sampler.m_data;
  };

  PlusMinusBasisFunctions basisFunctions{};
  basisFunctions.plusSide = getBasisFunctions(plusElementCoords);
  basisFunctions.minusSide = getBasisFunctions(minusElementCoords);

  return basisFunctions;
}

std::vector<double> getAllVertices(const seissol::dr::ReceiverPoints& receiverPoints) {
  std::vector<double> vertices(3 * (3 * receiverPoints.size()), 0.0);

  for (size_t pointIndex{0}; pointIndex < receiverPoints.size(); ++pointIndex) {
    for (int vertexIndex{0}; vertexIndex < ExtTriangle::size(); ++vertexIndex) {
      const auto& triangle = receiverPoints[pointIndex].globalTriangle;
      const auto& point = triangle.point(vertexIndex);

      const size_t globalVertexIndex = 3 * pointIndex + vertexIndex;
      for (int coordIndex{0}; coordIndex < ExtVrtxCoords::size(); ++coordIndex) {
        vertices[3 * globalVertexIndex + coordIndex] = point[coordIndex];
      }
    }
  }
  return vertices;
}

std::vector<unsigned int> getCellConnectivity(const seissol::dr::ReceiverPoints& receiverPoints) {
  std::vector<unsigned int> cells(3 * receiverPoints.size());

  for (size_t pointIndex{0}; pointIndex < receiverPoints.size(); ++pointIndex) {
    for (int vertexIndex{0}; vertexIndex < 3; ++vertexIndex) {
      const size_t globalVertexIndex = 3 * pointIndex + vertexIndex;
      cells[globalVertexIndex] = globalVertexIndex;
    }
  }
  return cells;
}

real computeTriangleArea(ExtTriangle& triangle) {
  const auto p0 = triangle.point(0).getAsEigen3LibVector();
  const auto p1 = triangle.point(1).getAsEigen3LibVector();
  const auto p2 = triangle.point(2).getAsEigen3LibVector();

  const auto vector1 = p1 - p0;
  const auto vector2 = p2 - p0;
  const auto normal = vector1.cross(vector2);
  return 0.5 * normal.norm();
}
} // namespace seissol::dr

namespace seissol::dr::filesystem_aux {
std::string getTimeStamp() {
  std::time_t time = std::time(nullptr);
  std::tm tm = *std::localtime(&time);

  std::stringstream timeStamp;
  timeStamp << std::put_time(&tm, "%F_%T");
  return timeStamp.str();
}

void generateBackupFileIfNecessary(std::string fileName, std::string fileExtension) {
  std::stringstream fullName;
  fullName << fileName << '.' << fileExtension;
  std::filesystem::path path(fullName.str());
  std::filesystem::directory_entry entry(path);

  if (entry.exists()) {
    auto stamp = getTimeStamp();
    std::stringstream backupFileName;
    backupFileName << fileName << ".bak_" << stamp << '.' << fileExtension;
    std::filesystem::path copyPath(backupFileName.str());
    std::filesystem::rename(path, copyPath);
  }
}
} // namespace seissol::dr::filesystem_aux
