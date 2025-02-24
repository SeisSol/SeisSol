// SPDX-FileCopyrightText: 2014-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#include "MeshTools.h"

#include <Geometry/MeshDefinition.h>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <vector>

namespace seissol {

const int MeshTools::FACE2NODES[4][3] = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
const int MeshTools::FACE2MISSINGNODE[4] = {3, 2, 1, 0};
const int MeshTools::NEIGHBORFACENODE2LOCAL[3] = {0, 2, 1};

void MeshTools::center(const Element& e, const std::vector<Vertex>& vertices, VrtxCoords center) {
  for (int i = 0; i < 3; i++) {
    center[i] = .25 * (vertices[e.vertices[0]].coords[i] + vertices[e.vertices[1]].coords[i] +
                       vertices[e.vertices[2]].coords[i] + vertices[e.vertices[3]].coords[i]);
  }
}

void MeshTools::center(const Element& e,
                       int face,
                       const std::vector<Vertex>& vertices,
                       VrtxCoords center) {
  for (int i = 0; i < 3; i++) {
    center[i] = (1.0 / 3.0) * (vertices[e.vertices[FACE2NODES[face][0]]].coords[i] +
                               vertices[e.vertices[FACE2NODES[face][1]]].coords[i] +
                               vertices[e.vertices[FACE2NODES[face][2]]].coords[i]);
  }
}

void MeshTools::normal(const Element& e,
                       int face,
                       const std::vector<Vertex>& vertices,
                       VrtxCoords normal) {
  VrtxCoords ab;
  VrtxCoords ac;
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
                                  VrtxCoords nrmal,
                                  VrtxCoords tangent1,
                                  VrtxCoords tangent2) {
  normal(e, face, vertices, nrmal);
  sub(vertices[e.vertices[FACE2NODES[face][1]]].coords,
      vertices[e.vertices[FACE2NODES[face][0]]].coords,
      tangent1);
  cross(nrmal, tangent1, tangent2);
}

void MeshTools::sub(const VrtxCoords v1, const VrtxCoords v2, VrtxCoords diff) {
  for (int i = 0; i < 3; i++) {
    diff[i] = v1[i] - v2[i];
  }
}

void MeshTools::mul(const VrtxCoords v, double s, VrtxCoords prod) {
  for (int i = 0; i < 3; i++) {
    prod[i] = v[i] * s;
  }
}

void MeshTools::cross(const VrtxCoords v1, const VrtxCoords v2, VrtxCoords cross) {
  cross[0] = v1[1] * v2[2] - v1[2] * v2[1];
  cross[1] = v1[2] * v2[0] - v1[0] * v2[2];
  cross[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

double MeshTools::dot(const VrtxCoords v1, const VrtxCoords v2) {
  double result = 0;
  for (int i = 0; i < 3; i++) {
    result += v1[i] * v2[i];
  }

  return result;
}

double MeshTools::norm(const VrtxCoords v) { return sqrt(norm2(v)); }

/**
 * Computes the square of the Euclidean norm
 */
double MeshTools::norm2(const VrtxCoords v) { return square(v[0]) + square(v[1]) + square(v[2]); }

double MeshTools::distance(const VrtxCoords v1, const VrtxCoords v2) {
  VrtxCoords diff;
  sub(v1, v2, diff);
  return norm(diff);
}

double MeshTools::surface(VrtxCoords faceNormal) {
  // Area of a triangle spanned by a and b = 0.5 * ||a x b||.
  return 0.5 * norm(faceNormal);
}

double MeshTools::surface(const Element& e, int face, const std::vector<Vertex>& vertices) {
  VrtxCoords n;
  normal(e, face, vertices, n);
  return surface(n);
}

double MeshTools::volume(const Element& e, const std::vector<Vertex>& vertices) {
  VrtxCoords ab;
  VrtxCoords ac;
  VrtxCoords ad;
  sub(vertices[e.vertices[1]].coords, vertices[e.vertices[0]].coords, ab);
  sub(vertices[e.vertices[2]].coords, vertices[e.vertices[0]].coords, ac);
  sub(vertices[e.vertices[3]].coords, vertices[e.vertices[0]].coords, ad);
  VrtxCoords area;
  cross(ab, ac, area);
  return fabs(dot(ad, area)) / 6.0;
}

void MeshTools::normalize(const VrtxCoords v, VrtxCoords vnormalized) {
  mul(v, 1.0 / norm(v), vnormalized);
}

void MeshTools::pointOnPlane(const Element& e,
                             int face,
                             const std::vector<Vertex>& vertices,
                             VrtxCoords result) {
  const size_t index = e.vertices[FACE2NODES[face][0]];
  assert(index < vertices.size());
  for (int i = 0; i < 3; i++) {
    result[i] = vertices[index].coords[i];
  }
}

bool MeshTools::inside(const Element& e, const std::vector<Vertex>& vertices, const VrtxCoords p) {
  VrtxCoords nrm;
  /* Our tetrahedron has 4 faces with the normals pointing outward.
   * The point is inside the tetrahedron if it lies on the backside
   * of each of the 4 planes defined by the normal vectors (and a point
   * on the plane). */
  for (unsigned face = 0; face < 4; ++face) {
    VrtxCoords pp;
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
