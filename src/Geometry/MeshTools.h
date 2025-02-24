// SPDX-FileCopyrightText: 2014-2024 SeisSol Group
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
  static void center(const Element& e, const std::vector<Vertex>& vertices, VrtxCoords center);

  /**
   * Computes the barycenter of a face of the element
   */
  static void
      center(const Element& e, int face, const std::vector<Vertex>& vertices, VrtxCoords center);

  /**
   * Computes the normal of a face
   */
  static void
      normal(const Element& e, int face, const std::vector<Vertex>& vertices, VrtxCoords normal);

  static void normalAndTangents(const Element& e,
                                int face,
                                const std::vector<Vertex>& vertices,
                                VrtxCoords nrmal,
                                VrtxCoords tangent1,
                                VrtxCoords tangent2);

  /**
   * Subtracts <code>v2</code> from <code>v1</code>
   */
  static void sub(const VrtxCoords v1, const VrtxCoords v2, VrtxCoords diff);

  /**
   * Multiplies <code>v</code> by the scalar <code>s</code>
   */
  static void mul(const VrtxCoords v, double s, VrtxCoords prod);

  /**
   * Computes the cross product of to vectors
   */
  static void cross(const VrtxCoords v1, const VrtxCoords v2, VrtxCoords cross);

  /**
   * Computes the dot product
   */
  static double dot(const VrtxCoords v1, const VrtxCoords v2);

  /**
   * Computes the Euclidean norm
   */
  static double norm(const VrtxCoords v);

  /**
   * Computes the square of the Euclidean norm
   */
  static double norm2(const VrtxCoords v);

  /**
   * Computes the Euclidean distance of two coordinates
   */
  static double distance(const VrtxCoords v1, const VrtxCoords v2);

  /**
   * Calculates the surface of a triangle based on its (unnormalized) normal.
   **/
  static double surface(VrtxCoords faceNormal);

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
  static void normalize(const VrtxCoords v, VrtxCoords vnormalized);

  /**
   * Returns a point on the plane spanned by the face-th plane.
   */
  static void pointOnPlane(const Element& e,
                           int face,
                           const std::vector<Vertex>& vertices,
                           VrtxCoords result);

  /**
   * Checks if a point p is inside a tetrahedron
   **/
  static bool inside(const Element& e, const std::vector<Vertex>& vertices, const VrtxCoords p);

  /** Maps from the face to the list of nodes */
  const static int FACE2NODES[4][3];

  /** Maps from the face to missing node of the element */
  const static int FACE2MISSINGNODE[4];

  /**
   * Maps from the neighbor face node id, to the local face node
   * Use <code>(3 + i - orientation) % 3</code> to respect the orientation
   * of the faces.
   */
  const static int NEIGHBORFACENODE2LOCAL[3];

  private:
  static double square(double v);
};

} // namespace seissol

#endif // SEISSOL_SRC_GEOMETRY_MESHTOOLS_H_
