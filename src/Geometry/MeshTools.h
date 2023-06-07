/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2014, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Useful mesh functions
 **/

#ifndef MESH_TOOLS_H
#define MESH_TOOLS_H

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

  static void normalAndTangents(Element const& e,
                                int face,
                                std::vector<Vertex> const& vertices,
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
  static double surface(Element const& e, int face, const std::vector<Vertex>& vertices);

  /**
   * Returns the volume of a tetrahedron.
   **/
  static double volume(Element const& e, const std::vector<Vertex>& vertices);

  /**
   * vnormalized = v / ||v||
   **/
  static void normalize(VrtxCoords const v, VrtxCoords vnormalized);

  /**
   * Returns a point on the plane spanned by the face-th plane.
   */
  static void pointOnPlane(Element const& e,
                           int face,
                           std::vector<Vertex> const& vertices,
                           VrtxCoords result);

  /**
   * Checks if a point p is inside a tetrahedron
   **/
  static bool inside(Element const& e, std::vector<Vertex> const& vertices, VrtxCoords const p);

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

#endif // MESH_TOOLS_H
