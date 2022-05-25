#include "OutputAux.hpp"
#include "Geometry/MeshTools.h"
#include "Numerical_aux/BasisFunction.h"
#include "Numerical_aux/Quadrature.h"
#include "Numerical_aux/Transformation.h"
#include <Eigen/Dense>
#include <limits>
#include <unordered_map>
#include <iomanip>
#include <filesystem>
#include <ctime>

namespace seissol::dr {

int getElementVertexId(int localSideId, int localFaceVertexId) {
  return MeshTools::FACE2NODES[localSideId][localFaceVertexId];
}

ExtTriangle getReferenceTriangle(int sideIdx) {
  ExtTriangle referenceFace;
  switch (sideIdx) {
  case 0:
    referenceFace.p1 = {0.0, 0.0, 0.0};
    referenceFace.p2 = {0.0, 1.0, 0.0};
    referenceFace.p3 = {1.0, 0.0, 0.0};
    break;
  case 1:
    referenceFace.p1 = {0.0, 0.0, 0.0};
    referenceFace.p2 = {1.0, 0.0, 0.0};
    referenceFace.p3 = {0.0, 0.0, 1.0};
    break;
  case 2:
    referenceFace.p1 = {0.0, 0.0, 0.0};
    referenceFace.p2 = {0.0, 0.0, 1.0};
    referenceFace.p3 = {0.0, 1.0, 0.0};
    break;
  case 3:
    referenceFace.p1 = {1.0, 0.0, 0.0};
    referenceFace.p2 = {0.0, 1.0, 0.0};
    referenceFace.p3 = {0.0, 0.0, 1.0};
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
    auto elementVertexId = getElementVertexId(localSideId, vertexId);
    auto globalVertexId = element.vertices[elementVertexId];

    triangle[vertexId] = verticesInfo[globalVertexId].coords;
  }
  return triangle;
}

ExtVrtxCoords getMidTrianglePoint(const ExtTriangle& triangle) {
  ExtVrtxCoords avgPoint{};
  for (int axis = 0; axis < 3; ++axis) {
    avgPoint.coords[axis] =
        (triangle.p1.coords[axis] + triangle.p2.coords[axis] + triangle.p3.coords[axis]) / 3.0;
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

std::tuple<unsigned, std::shared_ptr<double[]>, std::shared_ptr<double[]>>
    generateTriangleQuadrature(unsigned polyDegree) {

  // allocate data
  constexpr unsigned numQuadraturePoints = tensor::quadweights::Shape[0];
  std::shared_ptr<double[]> weights(new double[numQuadraturePoints],
                                    std::default_delete<double[]>());
  std::shared_ptr<double[]> points(new double[2 * numQuadraturePoints],
                                   std::default_delete<double[]>());

  // Generate triangle quadrature points and weights (Factory Method)
  auto pointsView = init::quadpoints::view::create(const_cast<real*>(init::quadpoints::Values));
  auto weightsView = init::quadweights::view::create(const_cast<real*>(init::quadweights::Values));
  auto* reshapedPoints = reshape<2>(&points[0]);
  for (size_t i = 0; i < numQuadraturePoints; ++i) {
    reshapedPoints[i][0] = pointsView(i, 0);
    reshapedPoints[i][1] = pointsView(i, 1);
    weights[i] = weightsView(i);
  }

  return std::make_tuple(numQuadraturePoints, weights, points);
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
    double nextPoint[2] = {facePoints[index][0], facePoints[index][1]};

    auto currentDistance = distance(targetPoint, nextPoint);
    if (shortestDistance > currentDistance) {
      shortestDistance = currentDistance;
      nearestPoint = static_cast<int>(index);
    }
  }
  return std::make_pair(nearestPoint, shortestDistance);
}

void assignNearestGaussianPoints(ReceiverPointsT& geoPoints) {
  std::shared_ptr<double[]> weights = nullptr;
  std::shared_ptr<double[]> pointsData = nullptr;
  unsigned numPoints{};

  std::tie(numPoints, weights, pointsData) = generateTriangleQuadrature(CONVERGENCE_ORDER + 1);
  double(*trianglePoints2D)[2] = reshape<2>(&pointsData[0]);

  for (auto& geoPoint : geoPoints) {

    double targetPoint2D[2];
    transformations::XiEtaZeta2chiTau(
        geoPoint.localFaceSideId, geoPoint.reference.coords, targetPoint2D);

    int nearestPoint{-1};
    double shortestDistance = std::numeric_limits<double>::max();
    std::tie(nearestPoint, shortestDistance) =
        getNearestFacePoint(targetPoint2D, trianglePoints2D, numPoints);
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
  auto distance = getDistanceFromPointToFace(point, face, faceNormal);
  double faceNormalLength = MeshTools::norm(faceNormal);
  auto adjustedDistance = distance / faceNormalLength;

  for (int i = 0; i < 3; ++i) {
    point.coords[i] += adjustedDistance * faceNormal[i];
  }
}

double getDistanceFromPointToFace(const ExtVrtxCoords& point,
                                  const ExtTriangle& face,
                                  const VrtxCoords faceNormal) {

  VrtxCoords diff{0.0, 0.0, 0.0};
  MeshTools::sub(face.p1.coords, point.coords, diff);

  // Note: faceNormal may not be precisely a unit vector
  double faceNormalLength = MeshTools::norm(faceNormal);
  return MeshTools::dot(faceNormal, diff) / faceNormalLength;
}

PlusMinusBasisFunctionsT getPlusMinusBasisFunctions(const VrtxCoords pointCoords,
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

  PlusMinusBasisFunctionsT basisFunctions{};
  basisFunctions.plusSide = getBasisFunctions(plusElementCoords);
  basisFunctions.minusSide = getBasisFunctions(minusElementCoords);

  return basisFunctions;
}

std::vector<double> getAllVertices(const seissol::dr::ReceiverPointsT& receiverPoints) {
  std::vector<double> vertices(3 * (3 * receiverPoints.size()), 0.0);

  for (size_t pointIndex{0}; pointIndex < receiverPoints.size(); ++pointIndex) {
    for (int vertexIndex{0}; vertexIndex < ExtTriangle::size(); ++vertexIndex) {
      const auto& triangle = receiverPoints[pointIndex].globalTriangle;
      auto& point = const_cast<ExtVrtxCoords&>(triangle.points[vertexIndex]);

      const size_t globalVertexIndex = 3 * pointIndex + vertexIndex;
      for (int coordIndex{0}; coordIndex < ExtVrtxCoords::size(); ++coordIndex) {
        vertices[3 * globalVertexIndex + coordIndex] = point[coordIndex];
      }
    }
  }
  return vertices;
}

std::vector<unsigned int> getCellConnectivity(const seissol::dr::ReceiverPointsT& receiverPoints) {
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
  auto p1 = triangle.p1.getAsEigenVector();
  auto p2 = triangle.p2.getAsEigenVector();
  auto p3 = triangle.p3.getAsEigenVector();

  auto vector1 = p2 - p1;
  auto vector2 = p3 - p1;
  auto normal = vector1.cross(vector2);
  return 0.5 * normal.norm();
}
} // namespace seissol::dr

namespace seissol::dr::os_support {
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
} // namespace seissol::dr::os_support
