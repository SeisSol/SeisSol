#ifndef SEISSOL_GEOMETRY_T_H
#define SEISSOL_GEOMETRY_T_H

#include "DynamicRupture/Output/OutputAux.hpp"
#include "DynamicRupture/Misc.h"
#include "Numerical_aux/BasisFunction.h"
#include "Numerical_aux/Transformation.h"
#include "Initializer/PointMapper.h"
#include "Geometry/MeshReader.h"
#include "tests/Geometry/MockReader.h"
#include <iostream>
#include <Eigen/Dense>
#include <tests/TestHelper.h>

namespace seissol::unit_test::dr {

using namespace seissol;
using namespace seissol::dr;

TEST_CASE("DR Geometry") {
  constexpr static int x{0};
  constexpr static int y{1};
  constexpr static int z{2};

  [[maybe_unused]] constexpr static int xi{0};
  constexpr static int eta{1};
  constexpr static int zeta{2};

  SUBCASE("Projection") {

    // Given a reference triangle in the first octant
    // TargetPoint - intersection of a line (which starts from the origin and goes along [1,1,1]
    // vector) and the inclined face (4th face)
    ExtVrtxCoords targetPoint{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};

    // 4th face
    ExtTriangle face(
        ExtVrtxCoords{1.0, 0.0, 0.0}, ExtVrtxCoords{0.0, 1.0, 0.0}, ExtVrtxCoords{0.0, 0.0, 1.0});

    constexpr double epsilon = 1e-6;
    {
      ExtVrtxCoords testPoint{0.0, 0.0, 0.0};
      VrtxCoords normalDirection{1.0, 1.0, 1.0};

      projectPointToFace(testPoint, face, normalDirection);

      REQUIRE(testPoint[x] == AbsApprox(targetPoint[x]).epsilon(epsilon));
      REQUIRE(testPoint[y] == AbsApprox(targetPoint[y]).epsilon(epsilon));
      REQUIRE(testPoint[z] == AbsApprox(targetPoint[z]).epsilon(epsilon));
    }
    {
      ExtVrtxCoords testPoint{1.0, 1.0, 1.0};
      VrtxCoords normalDirection{-1.0, -1.0, -1.0};

      projectPointToFace(testPoint, face, normalDirection);

      REQUIRE(testPoint[x] == AbsApprox(targetPoint[x]).epsilon(epsilon));
      REQUIRE(testPoint[y] == AbsApprox(targetPoint[y]).epsilon(epsilon));
      REQUIRE(testPoint[z] == AbsApprox(targetPoint[z]).epsilon(epsilon));
    }
  }

  SUBCASE("ClosestPoint") {
    double targetPoint[2] = {-0.25, -0.25};
    double facePoints[4][2] = {{1.0, 1.0}, {-1.0, 1.0}, {-1.0, -1.0}, {1.0, -1.0}};

    auto [testPointId, testDistance] = getNearestFacePoint(targetPoint, facePoints, 4);

    constexpr double epsilon = 1e-6;
    REQUIRE(testPointId == 2);
    REQUIRE(testDistance == AbsApprox(std::sqrt(2 * 0.75 * 0.75)).epsilon(epsilon));
  }

  SUBCASE("MiddlePoint") {
    ExtVrtxCoords point1{1.0, 2.0, 3.0};
    ExtVrtxCoords point2{-3.0, -2.0, -1.0};

    auto testMiddle = getMidPoint(point1, point2);

    constexpr double epsilon = 1e-6;
    REQUIRE(testMiddle[0] == AbsApprox(-1.0).epsilon(epsilon));
    REQUIRE(testMiddle[1] == AbsApprox(0.0).epsilon(epsilon));
    REQUIRE(testMiddle[2] == AbsApprox(1.0).epsilon(epsilon));
  }

  SUBCASE("MidTrianglePoint") {
    ExtVrtxCoords point1{0.5, 0.0, 2.0};
    ExtVrtxCoords point2{-0.5, 0.0, 3.0};
    ExtVrtxCoords point3{3.0, 1.0, -2.0};
    ExtTriangle triangle(point1, point2, point3);

    auto testMiddle = getMidPointTriangle(triangle);

    constexpr double epsilon = 1e-6;
    REQUIRE(testMiddle[x] == AbsApprox(1.0).epsilon(epsilon));
    REQUIRE(testMiddle[y] == AbsApprox(1 / 3.0).epsilon(epsilon));
    REQUIRE(testMiddle[z] == AbsApprox(1.0).epsilon(epsilon));
  }

  SUBCASE("TriangleQuadraturePoints") {
    // Coordinates are taken from the Fortran implementation
    double chiFortran[] = {
        0.94373743946307787,     0.94373743946307787,     0.94373743946307787,
        0.94373743946307787,     0.94373743946307787,     0.94373743946307787,
        0.94373743946307787,     0.81975930826310761,     0.81975930826310761,
        0.81975930826310761,     0.81975930826310761,     0.81975930826310761,
        0.81975930826310761,     0.81975930826310761,     0.64737528288683033,
        0.64737528288683033,     0.64737528288683033,     0.64737528288683033,
        0.64737528288683033,     0.64737528288683033,     0.64737528288683033,
        0.45284637366944464,     0.45284637366944464,     0.45284637366944464,
        0.45284637366944464,     0.45284637366944464,     0.45284637366944464,
        0.45284637366944464,     0.26578982278458946,     0.26578982278458946,
        0.26578982278458946,     0.26578982278458946,     0.26578982278458946,
        0.26578982278458946,     0.26578982278458946,     0.11467905316090421,
        0.11467905316090421,     0.11467905316090421,     0.11467905316090421,
        0.11467905316090421,     0.11467905316090421,     0.11467905316090421,
        2.2479386438712501E-002, 2.2479386438712501E-002, 2.2479386438712501E-002,
        2.2479386438712501E-002, 2.2479386438712501E-002, 2.2479386438712501E-002,
        2.2479386438712501E-002};

    double tauFortran[] = {
        5.4830900955589179E-002, 4.8991501878361855E-002, 3.9548223967454631E-002,
        2.8131280268461067E-002, 1.6714336569467501E-002, 7.2710586585602805E-003,
        1.4316595813329493E-003, 0.17565427919525450,     0.15694739278690259,
        0.12669525127960912,     9.0120345868446194E-002, 5.3545440457283260E-002,
        2.3293298949989799E-002, 4.5864125416378871E-003, 0.34365181310645293,
        0.30705347083287471,     0.24786787440468791,     0.17631235855658484,
        0.10475684270848173,     4.5571246280294943E-002, 8.9729040067167108E-003,
        0.53323073117395925,     0.47644255178423012,     0.38460663631768571,
        0.27357681316527771,     0.16254699001286968,     7.0711074546325303E-002,
        1.3922895156596097E-002, 0.71552743286656784,     0.63932496020254781,
        0.51609290886511228,     0.36710508860770530,     0.21811726835029835,
        9.4885217012862830E-002, 1.8682744348842751E-002, 0.86279303122343209,
        0.77090701909233450,     0.62231208026329454,     0.44266047341954790,
        0.26300886657580119,     0.11441392774676130,     2.2527915615663658E-002,
        0.95264658118522672,     0.85119131654161828,     0.68712130747329714,
        0.48876030678064375,     0.29039930608799031,     0.12632929701966925,
        2.4874032376060777E-002};

    auto data = generateTriangleQuadrature(7);
    double(*testTrianglePoints)[2] = unsafe_reshape<2>(&data.points[0]);

    constexpr double epsilon = 1e-6;
    for (unsigned i = 0; i < data.size; ++i) {
      REQUIRE(testTrianglePoints[i][0] == AbsApprox(chiFortran[i]).epsilon(epsilon));
      REQUIRE(testTrianglePoints[i][1] == AbsApprox(tauFortran[i]).epsilon(epsilon));
    }
  }

  SUBCASE("StrikeAndDipVectors") {
    VrtxCoords testNormal{-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
    VrtxCoords testStrike{0.0, 0.0, 0.0};
    VrtxCoords testDip{0.0, 0.0, 0.0};
    misc::computeStrikeAndDipVectors(testNormal, testStrike, testDip);

    // compute expected Strike results
    Eigen::Vector3d e3(0.0, 0.0, -1.0);
    Eigen::Vector3d normal(testNormal[0], testNormal[1], testNormal[2]);
    Eigen::Vector3d resultStrike = e3.cross(normal).normalized();

    constexpr double epsilon = 1e-6;
    for (unsigned i = 0; i < 3; ++i) {
      REQUIRE(testStrike[i] == AbsApprox(resultStrike(i)).epsilon(epsilon));
    }
    // compute expected Dip results
    Eigen::Vector3d resultDip = normal.cross(resultStrike);
    for (unsigned i = 0; i < 3; ++i) {
      REQUIRE(testDip[i] == AbsApprox(resultDip(i)).epsilon(epsilon));
    }
  }

  SUBCASE("XiEtaZeta2chiTau") {
    constexpr double epsilon = 1e-6;
    double testChiTau[2] = {0.0, 0.0};
    {
      unsigned face = 0;
      VrtxCoords xiEtaZeta{0.25, 0.1, 0.0};
      transformations::XiEtaZeta2chiTau(face, xiEtaZeta, testChiTau);
      REQUIRE(testChiTau[0] == AbsApprox(0.1).epsilon(epsilon));
      REQUIRE(testChiTau[1] == AbsApprox(0.25).epsilon(epsilon));
    }
    {
      unsigned face = 1;
      VrtxCoords xiEtaZeta{0.1, 0.0, 0.25};
      transformations::XiEtaZeta2chiTau(face, xiEtaZeta, testChiTau);
      REQUIRE(testChiTau[0] == AbsApprox(0.1).epsilon(epsilon));
      REQUIRE(testChiTau[1] == AbsApprox(0.25).epsilon(epsilon));
    }
    {
      unsigned face = 2;
      VrtxCoords xiEtaZeta{0.0, 0.1, 0.25};
      transformations::XiEtaZeta2chiTau(face, xiEtaZeta, testChiTau);
      REQUIRE(testChiTau[0] == AbsApprox(0.25).epsilon(epsilon));
      REQUIRE(testChiTau[1] == AbsApprox(0.1).epsilon(epsilon));
    }
    {
      unsigned face = 3;
      VrtxCoords xiEtaZeta{
          1 / 3.0, 1 / 3.0, 1 / 3.0}; // center of the 4th face (triangle in 3D space)
      transformations::XiEtaZeta2chiTau(face, xiEtaZeta, testChiTau);
      REQUIRE(testChiTau[0] == AbsApprox(1 / 3.0).epsilon(epsilon));
      REQUIRE(testChiTau[1] == AbsApprox(1 / 3.0).epsilon(epsilon));
    }
    {
      unsigned face = 3;
      ExtVrtxCoords xiEtaZeta{0.0, -0.15, 0.15};
      ExtTriangle fourthFace(
          ExtVrtxCoords{1.0, 0.0, 0.0}, ExtVrtxCoords{0.0, 1.0, 0.0}, ExtVrtxCoords{0.0, 0.0, 1.0});
      VrtxCoords normalDirection{1.0, 1.0, 1.0};
      projectPointToFace(xiEtaZeta, fourthFace, normalDirection);

      transformations::XiEtaZeta2chiTau(face, xiEtaZeta.coords, testChiTau);
      REQUIRE(testChiTau[0] == AbsApprox(xiEtaZeta[eta]).epsilon(epsilon));
      REQUIRE(testChiTau[1] == AbsApprox(xiEtaZeta[zeta]).epsilon(epsilon));
    }
  }

  SUBCASE("BasisFunctions") {

    VrtxCoords point{0.25, 0.25, 0.0};

    // placing two elements in such a way that basis functions on both sides end up being the same
    VrtxCoords plusElementCoords[4]{
        {2.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 2.0}};

    VrtxCoords minusElementCoords[4]{
        {2.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, -2.0}};

    const VrtxCoords* plusElementCoordsPtr[4]{
        &plusElementCoords[0], &plusElementCoords[1], &plusElementCoords[2], &plusElementCoords[3]};

    const VrtxCoords* minusElementCoordsPtr[4]{&minusElementCoords[0],
                                               &minusElementCoords[1],
                                               &minusElementCoords[2],
                                               &minusElementCoords[3]};

    auto basisFunctions =
        getPlusMinusBasisFunctions(point, plusElementCoordsPtr, minusElementCoordsPtr);

    constexpr double epsilon = 1e-6;
    for (unsigned i = 0; i < basisFunctions.plusSide.size(); ++i) {
      REQUIRE(basisFunctions.plusSide[i] ==
              AbsApprox(basisFunctions.minusSide[i]).epsilon(epsilon));
    }
  }

  SUBCASE("IsElementInside") {

    Eigen::Vector3d points[3] = {{0.25, 0.25, 0.25}, {0.5, 0.5, 0.5}, {0.75, 0.75, 0.1}};
    unsigned numPoints = 3;
    short contained[3] = {0, 0, 0};
    unsigned meshId[3] = {std::numeric_limits<unsigned>::max(),
                          std::numeric_limits<unsigned>::max(),
                          std::numeric_limits<unsigned>::max()};

    std::vector<Vertex> vertices;
    vertices.push_back({{0.0, 0.0, 0.0}, {0}});
    vertices.push_back({{1.0, 0.0, 0.0}, {0, 1}});
    vertices.push_back({{0.0, 1.0, 0.0}, {0, 1}});
    vertices.push_back({{0.0, 0.0, 1.0}, {0}});
    vertices.push_back({{1.0, 1.0, 0.0}, {1}});
    vertices.push_back({{1.0, 1.0, 1.0}, {1}});

    std::vector<Element> elements;
    Element e1;
    e1.localId = 0;
    e1.vertices[0] = 0;
    e1.vertices[1] = 1;
    e1.vertices[2] = 2;
    e1.vertices[3] = 3;
    elements.push_back(e1);

    Element e2;
    e2.localId = 1;
    e2.vertices[0] = 1;
    e2.vertices[1] = 4;
    e2.vertices[2] = 2;
    e2.vertices[3] = 5;
    elements.push_back(e2);

    initializers::findMeshIds(points, vertices, elements, numPoints, contained, meshId);

    REQUIRE(contained[0] == 1);
    REQUIRE(contained[1] == 0);
    REQUIRE(contained[2] == 1);

    REQUIRE(meshId[0] == 0);
    REQUIRE(meshId[1] == std::numeric_limits<unsigned>::max());
    REQUIRE(meshId[2] == 1);
  }
}
} // namespace seissol::unit_test::dr

#endif // SEISSOL_GEOMETRY_T_H