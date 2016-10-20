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
#ifndef GAMBITWRITER_H_
#define GAMBITWRITER_H_
#include <cstdio>
#include <ctime>

#include "GMSH.h"

template<unsigned DIM>
struct GambitInfo {
  static unsigned const Type;
};

template<unsigned DIM>
void writeNEU(char const* filename, GMSH<DIM> const& msh) {
  FILE* file;
  file = fopen(filename, "w");
  fprintf(file, "        CONTROL INFO 2.0.0\n");
  fprintf(file, "** GAMBIT NEUTRAL FILE\n");
  fprintf(file, "Gmsh mesh in GAMBIT neutral file format\n");
  fprintf(file, "PROGRAM:                Gambit     VERSION:  2.0.0\n");
  time_t rawtime;
  struct tm* timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  char datestring[256];
  strftime(datestring, sizeof(datestring), "%c", timeinfo);
  fprintf(file, "%s\n", datestring);
  fprintf(file, "     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL\n");
  fprintf(file, " %9d %9d %9d %9d %9d %9d\n", msh.numVertices, msh.numElements, msh.materialGroups.size(), msh.boundaryConditions.size(), DIM, DIM);
  fprintf(file, "ENDOFSECTION\n");
  
  // Vertices
  fprintf(file, "   NODAL COORDINATES 2.0.0\n");
  for (unsigned vertex = 0; vertex < msh.numVertices; ++vertex) {
    fprintf(file, "%10d", vertex+1);
    for (unsigned c = 0; c < DIM; ++c) {
      fprintf(file, "%20.11e", msh.vertices[vertex][c]);
    }
    fprintf(file, "\n");
  }
  fprintf(file, "ENDOFSECTION\n");
  
  // Elements
  fprintf(file, "      ELEMENTS/CELLS 2.0.0\n");
  for (unsigned element = 0; element < msh.numElements; ++element) {
    fprintf(file, "%8d %2d %2d ", element+1, GambitInfo<DIM>::Type, GMSH<DIM>::Element::NumNodes);
    for (unsigned n = 0; n < GMSH<DIM>::Element::NumNodes; ++n) {
      fprintf(file, "%8d", msh.elements[element].nodes[n]+1);
    }
    fprintf(file, "\n");
  }
  fprintf(file, "ENDOFSECTION\n");
  
  // Material groups
  for (auto group = msh.materialGroups.cbegin(); group != msh.materialGroups.cend(); ++group) {
    fprintf(file, "       ELEMENT GROUP 2.0.0\n");
    fprintf(file, "GROUP: %10d ELEMENTS: %10d MATERIAL: %10d NFLAGS: %10d\n", group->first, group->second.size(), 2, 1);
    fprintf(file, "Material group %d\n", group->first);
    fprintf(file, "       0");
    for (unsigned gm = 0; gm < group->second.size(); ++gm) {
      if (gm % 10 == 0) {
        fprintf(file, "\n");
      }
      fprintf(file, "%8d", group->second[gm]+1);
    }
    fprintf(file, "\n");
    fprintf(file, "ENDOFSECTION\n");
  }
  
  // Boundary conditions
  for (auto boundaryCondition = msh.boundaryConditions.cbegin(); boundaryCondition != msh.boundaryConditions.cend(); ++boundaryCondition) {
    fprintf(file, "       BOUNDARY CONDITIONS 2.0.0\n");
    // collect members in group
    fprintf(file, "%32d%8d%8d%8d%8d\n", boundaryCondition->first, 1, boundaryCondition->second.size(), 0, 6);
    for (auto boundary = boundaryCondition->second.begin(); boundary < boundaryCondition->second.end(); ++boundary) {
      fprintf(file, "%10d %5d %5d\n", boundary->element+1, GambitInfo<DIM>::Type, boundary->side+1);
    }
    fprintf(file, "ENDOFSECTION\n");
  }
  
  fclose(file);
}



#endif
