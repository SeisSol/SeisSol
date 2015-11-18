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

#include "GMSH.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>

#define MAX_NUMBER_OF_TAGS 16

struct Tri {
  unsigned nodes[3];
  unsigned boundaryCondition;
};

void error(std::string const& errMessage)
{
  std::cerr << "MSH parse error: " << errMessage << std::endl;
  exit(-1);
}

GMSH parseMSH(char const* filename)
{
  GMSH msh;
  msh.numVertices = 0;
  msh.numTetras = 0;

  Tri* tris = NULL;
  unsigned numTris = 0;
  
  std::ifstream in(filename, std::ifstream::in);
  
  unsigned tags[MAX_NUMBER_OF_TAGS];

  // Parse MSH
  while (in.good()) {
    std::string block;
    in >> block;
    if (block.compare("$Nodes") == 0) {
      in >> msh.numVertices;
      msh.vertices = new double[msh.numVertices][3];
      for (unsigned vertex = 0; vertex < msh.numVertices; ++vertex) {
        unsigned id;
        double x, y, z;
        in >> id >> x >> y >> z;
        if (id-1 != vertex) {
          error("Ids are not contiguous.");
        }
        msh.vertices[vertex][0] = x;
        msh.vertices[vertex][1] = y;
        msh.vertices[vertex][2] = z;
      }
    } else if (block.compare("$Elements") == 0) {
      unsigned numElements;
      in >> numElements;
      tris = new Tri[numElements];
      msh.tetras = new Tetra[numElements];
      
      for (unsigned element = 0; element < numElements; ++element) {
        unsigned id, type, numTags;
        in >> id >> type >> numTags;
        if (type != 2 && type != 4) {
          error("Unsuported element. Must be 2 (triangle) or 4 (tetrahedron).");
        }
        if (numTags > MAX_NUMBER_OF_TAGS) {
          error("Too many tags. Increase MAX_NUMBER_OF_TAGS.");
        }
        for (unsigned tag = 0; tag < numTags; ++tag) {
          in >> tags[tag];
        }
        if (type == 2) {
          Tri& tri = tris[numTris++];
          for (unsigned node = 0; node < 3; ++node) {
            in >> tri.nodes[node];
            --tri.nodes[node];
          }
          tri.boundaryCondition = tags[0];
        } else if (type == 4) {
          Tetra& tetra = msh.tetras[msh.numTetras];
          for (unsigned node = 0; node < 4; ++node) {
            in >> tetra.nodes[node];
            --tetra.nodes[node];
          }
          tetra.materialGroup = tags[0];
          msh.materialGroups[tetra.materialGroup].push_back(msh.numTetras);
          ++msh.numTetras;
        }
      }
    }
  }
  
  in.close();

  // vertex to triangle map
  std::vector<unsigned>* vertex2tri = vertex2tri = new std::vector<unsigned>[msh.numVertices];  
  for (unsigned tri = 0; tri < numTris; ++tri) {
    for (unsigned node = 0; node < 3; ++node) {
      vertex2tri[tris[tri].nodes[node]].push_back(tri);
    }
  }  
  for (unsigned vertex = 0; vertex < msh.numVertices; ++vertex) {
    std::sort(vertex2tri[vertex].begin(), vertex2tri[vertex].end());
  }
  
  // Uses GAMBIT neu conventions. See GAMBIT NEUTRAL FILE FORMAT Appendix C.2.
  unsigned const face2nodes[4][3] = {{1, 0, 2}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
  
  for (unsigned tetra = 0; tetra < msh.numTetras; ++tetra) {
    for (unsigned face = 0; face < 4; ++face) {
      unsigned nodes[3];
      for (unsigned node = 0; node < 3; ++node) {
        nodes[node] = msh.tetras[tetra].nodes[ face2nodes[face][node] ];
      }
      std::vector<unsigned> intersect01;
      std::vector<unsigned> intersect012;
      std::set_intersection( vertex2tri[ nodes[0] ].begin(), vertex2tri[ nodes[0] ].end(),
                             vertex2tri[ nodes[1] ].begin(), vertex2tri[ nodes[1] ].end(),
                             std::back_inserter(intersect01) );
      std::set_intersection( intersect01.begin(), intersect01.end(),
                             vertex2tri[ nodes[2] ].begin(), vertex2tri[ nodes[2] ].end(),
                             std::back_inserter(intersect012) );
      if (!intersect012.empty()) {
        if (intersect012.size() > 1) {
          error("A face of a tetrahedron exists multiple times in the surface mesh.");
        }
        Tri& tri = tris[ intersect012[0] ];
        Boundary bnd;
        bnd.element = tetra;
        bnd.side = face;
        msh.boundaryConditions[tri.boundaryCondition].push_back(bnd);
      }
    }
  }
  
  delete[] vertex2tri;
  delete[] tris;

  return msh;
}

