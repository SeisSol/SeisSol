/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
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
#include <cmath>

class MeshTools
{
public:
	/**
	 * Computes the barycenter of the element
	 */
	static void center(const Element& e, const std::vector<Vertex>& vertices, VrtxCoords center)
	{
		for (int i = 0; i < 3; i++) {
			center[i] = .25 * (vertices[e.vertices[0]].coords[i]
				+ vertices[e.vertices[1]].coords[i]
				+ vertices[e.vertices[2]].coords[i]
				+ vertices[e.vertices[3]].coords[i]);
		}
	}

	/**
	 * Computes the barycenter of a face of the element
	 */
	static void center(const Element& e, int face, const std::vector<Vertex>& vertices, VrtxCoords center)
	{
		for (int i = 0; i < 3; i++) {
			center[i] = (1.0/3.0) * (vertices[e.vertices[FACE2NODES[face][0]]].coords[i]
				+ vertices[e.vertices[FACE2NODES[face][1]]].coords[i]
				+ vertices[e.vertices[FACE2NODES[face][2]]].coords[i]);
		}
	}

	/**
	 * Computes the normal of a face
	 */
	static void normal(const Element& e, int face, const std::vector<Vertex>& vertices, VrtxCoords normal)
	{
		VrtxCoords ab;
		VrtxCoords ac;
		sub(vertices[e.vertices[FACE2NODES[face][1]]].coords, vertices[e.vertices[FACE2NODES[face][0]]].coords, ab);
		sub(vertices[e.vertices[FACE2NODES[face][2]]].coords, vertices[e.vertices[FACE2NODES[face][0]]].coords, ac);
		cross(ab, ac, normal);
	}
  
  static void normalAndTangents(Element const& e, int face, std::vector<Vertex> const& vertices, VrtxCoords nrmal, VrtxCoords tangent1, VrtxCoords tangent2)
  {
    normal(e, face, vertices, nrmal);
		sub(vertices[e.vertices[FACE2NODES[face][1]]].coords, vertices[e.vertices[FACE2NODES[face][0]]].coords, tangent1);
		cross(nrmal, tangent1, tangent2);
  }

	/**
	 * Subtracts <code>v2</code> from <code>v1</code>
	 */
	static void sub(const VrtxCoords v1, const VrtxCoords v2, VrtxCoords diff)
	{
		for (int i = 0; i < 3; i++)
			diff[i] = v1[i] - v2[i];
	}

	/**
	 * Multiplies <code>v</code> by the scalar <code>s</code>
	 */
	static void mul(const VrtxCoords v, double s, VrtxCoords prod)
	{
		for (int i = 0; i < 3; i++)
			prod[i] = v[i] * s;
	}

	/**
	 * Computes the cross product of to vectors
	 */
	static void cross(const VrtxCoords v1, const VrtxCoords v2, VrtxCoords cross)
	{
		cross[0] = v1[1]*v2[2] - v1[2]*v2[1];
		cross[1] = v1[2]*v2[0] - v1[0]*v2[2];
		cross[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}

	/**
	 * Computes the dot product
	 */
	static double dot(const VrtxCoords v1, const VrtxCoords v2)
	{
		double result = 0;
		for (int i = 0; i < 3; i++)
			result += v1[i]*v2[i];

		return result;
	}

	/**
	 * Computes the Euclidean norm
	 */
	static double norm(const VrtxCoords v)
	{
		return sqrt(norm2(v));
	}

	/**
	 * Computes the square of the Euclidean norm
	 */
	static double norm2(const VrtxCoords v)
	{
		return square(v[0]) + square(v[1]) + square(v[2]);
	}

	/**
	 * Computes the Euclidean distance of two coordinates
	 */
	static double distance(const VrtxCoords v1, const VrtxCoords v2)
	{
		VrtxCoords diff;
		sub(v1, v2, diff);
		return norm(diff);
	}
  
  /**
   * Calculates the surface of a triangle based on its (unnormalized) normal.
   **/
  static double surface(VrtxCoords faceNormal)
  {
    // Area of a triangle spanned by a and b = 0.5 * ||a x b||.
    return 0.5 * norm(faceNormal);
  }
  
  /**
   * Returns the surface area of the side of a tetrahedron.
   **/
  static double surface(Element const& e, int face, const std::vector<Vertex>& vertices)
  {
    VrtxCoords N;
    normal(e, face, vertices, N);
    return surface(N);
  }
  
  /**
   * Returns the volume of a tetrahedron.
   **/
  static double volume(Element const& e, const std::vector<Vertex>& vertices)
  {
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
  
  /**
   * vnormalized = v / ||v||
   **/
	static void normalize(VrtxCoords const v, VrtxCoords vnormalized)
  {
    mul(v, 1.0 / norm(v), vnormalized);
  }
  
  /**
   * Returns a point on the plane spanned by the face-th plane.
   */
  static void pointOnPlane(Element const& e, int face, std::vector<Vertex> const& vertices, VrtxCoords result)
  {
    int index = e.vertices[FACE2NODES[face][0]];
    for (int i = 0; i < 3; i++) {
			result[i] = vertices[index].coords[i];
    }
  }
  
  /**
   * Checks if a point p is inside a tetrahedron
   **/
  static bool inside(Element const&e, std::vector<Vertex> const& vertices, VrtxCoords const p)
  {
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
	static double square(double v)
	{
		return v*v;
	}
};

#endif // MESH_TOOLS_H
