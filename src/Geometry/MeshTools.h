// SPDX-FileCopyrightText: 2014 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_GEOMETRY_MESHTOOLS_H_
#define SEISSOL_SRC_GEOMETRY_MESHTOOLS_H_

#include "MeshDefinition.h"

#include <vector>

namespace seissol {

class MeshTools {
  public:
  /**
   * Computes the barycenter of the element
   */
  static void center(const Element& e, const std::vector<Vertex>& vertices, CoordinateT& center);

  /**
   * Computes the barycenter of a face of the element
   */
  static void
      center(const Element& e, int face, const std::vector<Vertex>& vertices, CoordinateT& center);

  /**
   * Computes the normal of a face
   */
  static void
      normal(const Element& e, int face, const std::vector<Vertex>& vertices, CoordinateT& normal);

  static void normalAndTangents(const Element& e,
                                int face,
                                const std::vector<Vertex>& vertices,
                                CoordinateT& nrmal,
                                CoordinateT& tangent1,
                                CoordinateT& tangent2);

  /**
   * Subtracts <code>v2</code> from <code>v1</code>
   */
  static void sub(const CoordinateT& v1, const CoordinateT& v2, CoordinateT& diff);

  /**
   * Multiplies <code>v</code> by the scalar <code>s</code>
   */
  static void mul(const CoordinateT& v, double s, CoordinateT& prod);

  /**
   * Computes the cross product of to vectors
   */
  static void cross(const CoordinateT& v1, const CoordinateT& v2, CoordinateT& cross);

  /**
   * Computes the dot product
   */
  static double dot(const CoordinateT& v1, const CoordinateT& v2);

  /**
   * Computes the Euclidean norm
   */
  static double norm(const CoordinateT& v);

  /**
   * Computes the square of the Euclidean norm
   */
  static double norm2(const CoordinateT& v);

  /**
   * Computes the Euclidean distance of two coordinates
   */
  static double distance(const CoordinateT& v1, const CoordinateT& v2);

  /**
   * Calculates the surface of a triangle based on its (unnormalized) normal.
   **/
  static double surface(CoordinateT& faceNormal);

  /**
   * Returns the surface area of the side of a tetrahedron.
   **/
  static double surface(const Element& e, int face, const std::vector<Vertex>& vertices);

  /**
   * Returns the volume of a tetrahedron.
   **/
  static double volume(const Element& e, const std::vector<Vertex>& vertices);

  /**
   * vnormalized = v / ||v||
   **/
  static void normalize(const CoordinateT& v, CoordinateT& vnormalized);

  /**
   * Returns a point on the plane spanned by the face-th plane.
   */
  static void pointOnPlane(const Element& e,
                           int face,
                           const std::vector<Vertex>& vertices,
                           CoordinateT& result);

  /**
   * Checks if a point p is inside a tetrahedron
   **/
  static bool inside(const Element& e, const std::vector<Vertex>& vertices, const CoordinateT& p);

  /** Maps from the face to the list of nodes */
  const static std::array<std::array<int, 3>, 4> FACE2NODES;

  /** Maps from the face to missing node of the element */
  const static std::array<int, 4> FACE2MISSINGNODE;

  /**
   * Maps from the neighbor face node id, to the local face node
   * Use <code>(3 + i - orientation) % 3</code> to respect the orientation
   * of the faces.
   */
  const static std::array<int, 3> NEIGHBORFACENODE2LOCAL;

  private:
  static double square(double v);
};

} // namespace seissol

#endif // SEISSOL_SRC_GEOMETRY_MESHTOOLS_H_
