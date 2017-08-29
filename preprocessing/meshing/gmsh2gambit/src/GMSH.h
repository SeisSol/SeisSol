/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
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
 **/
#pragma once

#include <fstream>
#include <algorithm>
#include <vector>
#include <unordered_map>

namespace
{
  #define MAX_NUMBER_OF_TAGS 16
}

void error(std::string const& errMessage);
unsigned boundary_code(std::string name);

template<unsigned N>
struct Simplex {
  /* According to MSH mesh format documentation:
   * 1 = 2-node line.
   * 2 = 3-node triangle.
   * 4 = 4-node tetrahedron.
   * Will not work for other types (e.g. 11 = 10-node second order tetrahedron).
   */
  static unsigned const Type = N <= 1 ? 1 : 1 << (N-1);
  static unsigned const NumNodes = N+1;
  static unsigned const NumFaces = N+1;
  /* Each N-simplex includes N+1 (N-1)-simplices where each has DIM nodes
   * E.g. Triangle    (N=2): 3 lines with 2 nodes each
   *      Tetrahedron (N=3): 4 triangles with 3 nodes each */
  static unsigned const Face2Nodes[N+1][N];

  unsigned nodes[NumNodes];
  unsigned group;
};

struct Boundary {
  unsigned element;
  unsigned side;
};

struct Region {
  explicit Region(unsigned id, unsigned type, std::string name);
  explicit Region();
  virtual ~Region() = default;

  unsigned id;
  unsigned type;
  std::string name;
};

template<unsigned DIM>
struct GMSH {
  typedef Simplex<DIM> Tetrahedron;
  typedef Simplex<DIM-1> Face;
  typedef Simplex<DIM-2> Line;

  double (*vertices)[3];
  Tetrahedron* tetrahedra;

  unsigned numVertices;
  unsigned numTetrahedra;

  std::unordered_map<unsigned, std::vector<unsigned> > materialGroups;
  std::unordered_map<unsigned, std::vector<Boundary> > boundaryConditions;
  std::unordered_map<unsigned, Region> regions;

  explicit GMSH(char const* filename);

  virtual ~GMSH() {
    delete[] vertices;
    delete[] tetrahedra;
  }
};

template<unsigned DIM>
GMSH<DIM>::GMSH(char const* filename)
{
  numVertices = 0;
  numTetrahedra = 0;

  Face* faces = NULL;
  unsigned numFaces = 0;

  std::ifstream in(filename, std::ifstream::in);

  unsigned tags[MAX_NUMBER_OF_TAGS];

  // Parse MSH
  while (in.good()) {
    std::string block;
    in >> block;
    if (block.compare("$PhysicalNames") == 0) {
      unsigned numRegions;
      in >> numRegions;

      for (auto i = 0; i < numRegions; ++i) {
        unsigned dim, id;
        std::string name;
        in >> dim >> id >> name;
        name = name.substr(1, name.size()-2);

        regions[id] = Region(id, boundary_code(name), name);
      }
    } else if (block.compare("$Nodes") == 0) {
      in >> numVertices;
      vertices = new double[numVertices][3];
      for (unsigned vertex = 0; vertex < numVertices; ++vertex) {
        unsigned id;
        double x, y, z;
        in >> id >> x >> y >> z;
        if (id-1 != vertex) {
          error("Ids are not contiguous.");
        }
        vertices[vertex][0] = x;
        vertices[vertex][1] = y;
        vertices[vertex][2] = z;
      }
    } else if (block.compare("$Elements") == 0) {
      unsigned numGeoShapes;
      in >> numGeoShapes;
      faces = new Face[numGeoShapes];
      tetrahedra = new Tetrahedron[numGeoShapes];

      for (unsigned geoShape = 0; geoShape < numGeoShapes; ++geoShape) {
        unsigned id, type, numTags;
        in >> id >> type >> numTags;
        if (numTags > MAX_NUMBER_OF_TAGS) {
          error("Too many tags. Increase MAX_NUMBER_OF_TAGS.");
        }
        for (unsigned tag = 0; tag < numTags; ++tag) {
          in >> tags[tag];
        }
        if (type == Tetrahedron::Type) {
          Tetrahedron& tet = tetrahedra[numTetrahedra];
          for (unsigned node = 0; node < Tetrahedron::NumNodes; ++node) {
            in >> tet.nodes[node];
            --tet.nodes[node];
          }
          tet.group = tags[0];
          materialGroups[tet.group].push_back(numTetrahedra);
          ++numTetrahedra;
        } else if (type == Face::Type) {
          Face& face = faces[numFaces++];
          for (unsigned node = 0; node < Face::NumNodes; ++node) {
            in >> face.nodes[node];
            --face.nodes[node];
          }
          face.group = tags[0];
        } else if (type == Line::Type) {
          unsigned dummy;
          for (unsigned node = 0; node < Line::NumNodes; ++node) {
            in >> dummy;
          }
        } else {
          error("Unsupported type. Must be 1 (line, 2D), 2 (triangle, 2D and 3D), or 4 (tetrahedron, 3D).");
        }
      }
    }
  }

  in.close();

  // vertex to triangle map
  std::vector<unsigned>* vertex2face = vertex2face = new std::vector<unsigned>[numVertices];
  for (unsigned face = 0; face < numFaces; ++face) {
    for (unsigned node = 0; node < Face::NumNodes; ++node) {
      vertex2face[faces[face].nodes[node]].push_back(face);
    }
  }
  for (unsigned vertex = 0; vertex < numVertices; ++vertex) {
    std::sort(vertex2face[vertex].begin(), vertex2face[vertex].end());
  }

  for (unsigned tet = 0; tet < numTetrahedra; ++tet) {
    for (unsigned elmtFace = 0; elmtFace < Tetrahedron::NumFaces; ++elmtFace) {
      unsigned nodes[Face::NumNodes];
      for (unsigned node = 0; node < Face::NumNodes; ++node) {
        nodes[node] = tetrahedra[tet].nodes[Tetrahedron::Face2Nodes[elmtFace][node]];
      }
      std::vector<unsigned> intersect[Face::NumNodes-1];
      std::set_intersection( vertex2face[ nodes[0] ].begin(), vertex2face[ nodes[0] ].end(),
                             vertex2face[ nodes[1] ].begin(), vertex2face[ nodes[1] ].end(),
                             std::back_inserter(intersect[0]) );
      for (unsigned node = 2; node < Face::NumNodes; ++node) {
        std::set_intersection( intersect[node-2].begin(), intersect[node-2].end(),
                               vertex2face[ nodes[node] ].begin(), vertex2face[ nodes[node] ].end(),
                             std::back_inserter(intersect[node-1]) );
      }
      if (!intersect[Face::NumNodes-2].empty()) {
        if (intersect[Face::NumNodes-2].size() > 1) {
          error("A face of an tet exists multiple times in the surface mesh.");
        }
        Face& face = faces[ intersect[Face::NumNodes-2][0] ];
        Boundary bnd;
        bnd.element = tet;
        bnd.side = elmtFace;
        boundaryConditions[face.group].push_back(bnd);
      }
    }
  }

  delete[] vertex2face;
  delete[] faces;
}
