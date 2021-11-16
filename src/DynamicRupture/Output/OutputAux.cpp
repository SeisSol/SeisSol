#include "OutputAux.hpp"
#include "DynamicRupture/Math.h"
#include "Geometry/MeshTools.h"
#include "Numerical_aux/Quadrature.h"
#include "Numerical_aux/Transformation.h"
#include "Numerical_aux/BasisFunction.h"
#include <Eigen/Dense>
#include <limits>

namespace seissol::dr {

    int getElementVertexId(int localSideId, int localFaceVertexId) {
      return MeshTools::FACE2NODES[localSideId][localFaceVertexId];
    }


    ExtTriangle getReferenceFace(int localSideId) {
      ExtTriangle referenceFace;
      switch (localSideId) {
        case 0:
          referenceFace.p1.xi = 0.0; referenceFace.p1.eta = 0.0; referenceFace.p1.zeta = 0.0;
          referenceFace.p2.xi = 0.0; referenceFace.p2.eta = 1.0; referenceFace.p2.zeta = 0.0;
          referenceFace.p3.xi = 1.0; referenceFace.p3.eta = 0.0; referenceFace.p3.zeta = 0.0;
          break;
        case 1:
          referenceFace.p1.xi = 0.0; referenceFace.p1.eta = 0.0; referenceFace.p1.zeta = 0.0;
          referenceFace.p2.xi = 1.0; referenceFace.p2.eta = 0.0; referenceFace.p2.zeta = 0.0;
          referenceFace.p3.xi = 0.0; referenceFace.p3.eta = 0.0; referenceFace.p3.zeta = 1.0;
          break;
        case 2:
          referenceFace.p1.xi = 0.0; referenceFace.p1.eta = 0.0; referenceFace.p1.zeta = 0.0;
          referenceFace.p2.xi = 0.0; referenceFace.p2.eta = 0.0; referenceFace.p2.zeta = 1.0;
          referenceFace.p3.xi = 0.0; referenceFace.p3.eta = 1.0; referenceFace.p3.zeta = 0.0;
          break;
        case 3:
          referenceFace.p1.xi = 1.0; referenceFace.p1.eta = 0.0; referenceFace.p1.zeta = 0.0;
          referenceFace.p2.xi = 0.0; referenceFace.p2.eta = 1.0; referenceFace.p2.zeta = 0.0;
          referenceFace.p3.xi = 0.0; referenceFace.p3.eta = 0.0; referenceFace.p3.zeta = 1.0;
          break;
        default:
          throw std::runtime_error("Unknown Local Side Id. Must be 0, 1, 2 or 3");
      }

      return referenceFace;
    }


    void computeStrikeAndDipVectors(const VrtxCoords normal, VrtxCoords strike, VrtxCoords dip) {
      // Note: equations are explained in documentation -> left-lateral-right-lateral-normal-reverse

      // compute normalized strike vector
      auto StrikeInvLength = 1.0 / std::sqrt(normal[0] * normal[0] + normal[1] * normal[1]);
      strike[0] = normal[1] * StrikeInvLength;
      strike[1] = -normal[0] * StrikeInvLength;
      strike[2] = 0.0;

      // compute normalized dip vector
      dip[0] = -1.0 * strike[1] * normal[2];
      dip[1] = strike[0] * normal[2];
      dip[2] = strike[1] * normal[0] - strike[0] * normal[1];
      auto dipInvLength = 1.0 / std::sqrt(dip[0] * dip[0] + dip[1] * dip[1] + dip[2] * dip[2]);
      dip[0] *= dipInvLength;
      dip[1] *= dipInvLength;
      dip[2] *= dipInvLength;
    }


    ExtVrtxCoords getMidTrianglePoint(const ExtTriangle& triangle) {
      ExtVrtxCoords avgPoint{};
      for (int axis = 0; axis < 3; ++axis) {
        avgPoint.coords[axis] = (triangle.p1.coords[axis] + triangle.p2.coords[axis] + triangle.p3.coords[axis]) / 3.0;
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


    std::tuple<unsigned, std::shared_ptr<double []>, std::shared_ptr<double []>>
    generateTriangleQuadrature(unsigned polyDegree) {

      // allocate data
      unsigned numQuadraturePoints = polyDegree * polyDegree;
      std::shared_ptr<double []> weights( new double[numQuadraturePoints], std::default_delete<double[]>());
      std::shared_ptr<double []> points( new double[2 * numQuadraturePoints], std::default_delete<double[]>());

      // Generate triangle quadrature points and weights (Factory Method)
      seissol::quadrature::TriangleQuadrature(reshape<2>(&points[0]), &weights[0], polyDegree);
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
          nearestPoint = index;
        }
      }
      return std::make_pair(nearestPoint, shortestDistance);
    }


    void assignNearestGaussianPoints(ReceiverPointsT& geoPoints) {
      std::shared_ptr<double []> weights = nullptr;
      std::shared_ptr<double []> pointsData = nullptr;
      unsigned numPoints{};

      std::tie(numPoints, weights, pointsData) = generateTriangleQuadrature(CONVERGENCE_ORDER + 1);
      double (*TrianglePoints2D)[2] = reshape<2>(&pointsData[0]);

      for (auto& geoPoint: geoPoints) {

        double targetPoint2D[2];
        transformations::XiEtaZeta2chiTau(geoPoint.localFaceSideId, geoPoint.referece.coords, targetPoint2D);

        int nearestPoint{-1};
        double shortestDistance = std::numeric_limits<double>::max();
        std::tie(nearestPoint, shortestDistance) = getNearestFacePoint(targetPoint2D,
                                                                       TrianglePoints2D,
                                                                       numPoints);
        geoPoint.nearestGpIndex = nearestPoint;
        geoPoint.distanceToNearestGp = shortestDistance;
      }
    }


    void projectPointToFace(ExtVrtxCoords& point, const ExtTriangle& face, const VrtxCoords faceNormal) {
      using namespace Eigen;

      Vector3d originalPoint(point.x, point.y, point.z);

      Vector3d r = originalPoint - Vector3d(face.p1.x, face.p1.y, face.p1.z);
      Vector3d direction(faceNormal[0], faceNormal[1], faceNormal[2]);
      direction.normalize();

      auto displacement = -direction.dot(r);
      Vector3d targetPoint = originalPoint + displacement * direction;

      for (int i = 0; i < 3; ++i) {
        point.coords[i] = targetPoint(i);
      }
    }

    PlusMinusBasisFunctionsT getPlusMinusBasisFunctions(const VrtxCoords pointCoords,
                                                        const VrtxCoords* plusElementCoords[4],
                                                        const VrtxCoords* minusElementCoords[4]) {

      PlusMinusBasisFunctionsT basisFunctions{};
      Eigen::Vector3d point(pointCoords[0], pointCoords[1], pointCoords[2]);

      {
        auto plusXiEtaZeta = transformations::tetrahedronGlobalToReference(*plusElementCoords[0],
                                                                           *plusElementCoords[1],
                                                                           *plusElementCoords[2],
                                                                           *plusElementCoords[3],
                                                                           point);

        basisFunction::SampledBasisFunctions<real> sampler(CONVERGENCE_ORDER,
                                                           plusXiEtaZeta[0],
                                                           plusXiEtaZeta[1],
                                                           plusXiEtaZeta[2]);
        basisFunctions.plusSide = std::move(sampler.m_data);
      }

      {
        auto minusXiEtaZeta = transformations::tetrahedronGlobalToReference(*minusElementCoords[0],
                                                                            *minusElementCoords[1],
                                                                            *minusElementCoords[2],
                                                                            *minusElementCoords[3],
                                                                            point);

        basisFunction::SampledBasisFunctions<real> sampler(CONVERGENCE_ORDER,
                                                           minusXiEtaZeta[0],
                                                           minusXiEtaZeta[1],
                                                           minusXiEtaZeta[2]);
        basisFunctions.minusSide = std::move(sampler.m_data);
      }

      return basisFunctions;
    }


    std::vector<double> getAllVertices(const seissol::dr::ReceiverPointsT& receiverPoints) {
      std::vector<double> vertices(3 * (3  * receiverPoints.size()), 0.0);

      for (size_t pointIndex{0}; pointIndex < receiverPoints.size(); ++pointIndex) {
        for (int vertexIndex{0}; vertexIndex < 3; ++vertexIndex) {
          const size_t globalVertexIndex = 3 * pointIndex + vertexIndex;

          vertices[3 * globalVertexIndex] = receiverPoints[pointIndex].globalSubTet.points[vertexIndex].x;
          vertices[3 * globalVertexIndex + 1] = receiverPoints[pointIndex].globalSubTet.points[vertexIndex].y;
          vertices[3 * globalVertexIndex + 2] = receiverPoints[pointIndex].globalSubTet.points[vertexIndex].z;
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
}