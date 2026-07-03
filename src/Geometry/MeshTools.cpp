// SPDX-FileCopyrightText: 2014 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#include "MeshTools.h"

#include "Common/Constants.h"
#include "Geometry/MeshDefinition.h"

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <vector>

namespace seissol {

const std::array<std::array<int, 3>, 4> MeshTools::FACE2NODES = {std::array<int, 3>{0, 2, 1},
                                                                 std::array<int, 3>{0, 1, 3},
                                                                 std::array<int, 3>{0, 3, 2},
                                                                 std::array<int, 3>{1, 2, 3}};
const std::array<int, 4> MeshTools::FACE2MISSINGNODE = {3, 2, 1, 0};
const std::array<int, 3> MeshTools::NEIGHBORFACENODE2LOCAL = {0, 2, 1};

void MeshTools::center(const Element& e, const std::vector<Vertex>& vertices, CoordinateT& center) {
  for (std::size_t i = 0; i < Cell::Dim; i++) {
    center[i] = .25 * (vertices[e.vertices[0]].coords[i] + vertices[e.vertices[1]].coords[i] +
                       vertices[e.vertices[2]].coords[i] + vertices[e.vertices[3]].coords[i]);
  }
}

void MeshTools::center(const Element& e,
                       int face,
                       const std::vector<Vertex>& vertices,
                       CoordinateT& center) {
  for (std::size_t i = 0; i < Cell::Dim; i++) {
    center[i] = (1.0 / 3.0) * (vertices[e.vertices[FACE2NODES[face][0]]].coords[i] +
                               vertices[e.vertices[FACE2NODES[face][1]]].coords[i] +
                               vertices[e.vertices[FACE2NODES[face][2]]].coords[i]);
  }
}

void MeshTools::normal(const Element& e,
                       int face,
                       const std::vector<Vertex>& vertices,
                       CoordinateT& normal) {
  CoordinateT ab;
  CoordinateT ac;
  sub(vertices[e.vertices[FACE2NODES[face][1]]].coords,
      vertices[e.vertices[FACE2NODES[face][0]]].coords,
      ab);
  sub(vertices[e.vertices[FACE2NODES[face][2]]].coords,
      vertices[e.vertices[FACE2NODES[face][0]]].coords,
      ac);
  cross(ab, ac, normal);
}

void MeshTools::normalAndTangents(const Element& e,
                                  int face,
                                  const std::vector<Vertex>& vertices,
                                  CoordinateT& nrmal,
                                  CoordinateT& tangent1,
                                  CoordinateT& tangent2) {
  normal(e, face, vertices, nrmal);
  sub(vertices[e.vertices[FACE2NODES[face][1]]].coords,
      vertices[e.vertices[FACE2NODES[face][0]]].coords,
      tangent1);
  cross(nrmal, tangent1, tangent2);
}

void MeshTools::sub(const CoordinateT& v1, const CoordinateT& v2, CoordinateT& diff) {
  for (std::size_t i = 0; i < Cell::Dim; i++) {
    diff[i] = v1[i] - v2[i];
  }
}

void MeshTools::mul(const CoordinateT& v, double s, CoordinateT& prod) {
  for (std::size_t i = 0; i < Cell::Dim; i++) {
    prod[i] = v[i] * s;
  }
}

void MeshTools::cross(const CoordinateT& v1, const CoordinateT& v2, CoordinateT& cross) {
  cross[0] = v1[1] * v2[2] - v1[2] * v2[1];
  cross[1] = v1[2] * v2[0] - v1[0] * v2[2];
  cross[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

double MeshTools::dot(const CoordinateT& v1, const CoordinateT& v2) {
  double result = 0;
  for (int i = 0; i < 3; i++) {
    result += v1[i] * v2[i];
  }

  return result;
}

double MeshTools::norm(const CoordinateT& v) { return sqrt(norm2(v)); }

/**
 * Computes the square of the Euclidean norm
 */
double MeshTools::norm2(const CoordinateT& v) { return square(v[0]) + square(v[1]) + square(v[2]); }

double MeshTools::distance(const CoordinateT& v1, const CoordinateT& v2) {
  CoordinateT diff;
  sub(v1, v2, diff);
  return norm(diff);
}

double MeshTools::surface(CoordinateT& faceNormal) {
  // Area of a triangle spanned by a and b = 0.5 * ||a x b||.
  return 0.5 * norm(faceNormal);
}

double MeshTools::surface(const Element& e, int face, const std::vector<Vertex>& vertices) {
  CoordinateT n;
  normal(e, face, vertices, n);
  return surface(n);
}

double MeshTools::volume(const Element& e, const std::vector<Vertex>& vertices) {
  CoordinateT ab;
  CoordinateT ac;
  CoordinateT ad;
  sub(vertices[e.vertices[1]].coords, vertices[e.vertices[0]].coords, ab);
  sub(vertices[e.vertices[2]].coords, vertices[e.vertices[0]].coords, ac);
  sub(vertices[e.vertices[3]].coords, vertices[e.vertices[0]].coords, ad);
  CoordinateT area;
  cross(ab, ac, area);
  return fabs(dot(ad, area)) / 6.0;
}

void MeshTools::normalize(const CoordinateT& v, CoordinateT& vnormalized) {
  mul(v, 1.0 / norm(v), vnormalized);
}

void MeshTools::pointOnPlane(const Element& e,
                             int face,
                             const std::vector<Vertex>& vertices,
                             CoordinateT& result) {
  const size_t index = e.vertices[FACE2NODES[face][0]];
  assert(index < vertices.size());
  for (std::size_t i = 0; i < Cell::Dim; i++) {
    result[i] = vertices[index].coords[i];
  }
}

bool MeshTools::inside(const Element& e,
                       const std::vector<Vertex>& vertices,
                       const CoordinateT& p) {
  CoordinateT nrm;
  /* Our tetrahedron has 4 faces with the normals pointing outward.
   * The point is inside the tetrahedron if it lies on the backside
   * of each of the 4 planes defined by the normal vectors (and a point
   * on the plane). */
  static_assert(Cell::NumFaces == 4, "Non-tetrahedral meshes are not supported here yet.");
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    CoordinateT pp;
    sub(p, vertices[e.vertices[FACE2NODES[face][0]]].coords, pp);
    normal(e, face, vertices, nrm);
    if (dot(nrm, pp) > 0.0) {
      return false;
    }
  }

  return true;
}

double MeshTools::square(double v) { return v * v; }

} // namespace seissol
