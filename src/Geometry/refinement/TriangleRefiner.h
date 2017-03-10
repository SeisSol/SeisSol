/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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
 */

#ifndef TRIANGLE_REFINER_H
#define TRIANGLE_REFINER_H

#include <vector>
#include <glm/vec2.hpp>
#include <cmath>

namespace seissol {
  namespace refinement {
    struct Triangle;
    class TriangleRefiner;
  }
}

struct seissol::refinement::Triangle
{
  glm::dvec2 x[3];
  double area;
};

class seissol::refinement::TriangleRefiner
{
private:  
  void refine(glm::dvec2 x0, glm::dvec2 x1, glm::dvec2 x2, unsigned depth, unsigned maxDepth)
  {
    glm::dvec2 chi = x1-x0;
    glm::dvec2 tau = x2-x0;
      
    if (depth < maxDepth) {
      // Edge midpoints
      glm::dvec2 m0 = 0.5 * chi + 0.5 * tau + x0;
      glm::dvec2 m1 = 0.5 * tau + x0;
      glm::dvec2 m2 = 0.5 * chi + x0;
      
      refine(x0, m2, m1, depth+1, maxDepth);
      refine(m2, x1, m0, depth+1, maxDepth);
      refine(m1, m0, x2, depth+1, maxDepth);
      refine(m0, m1, m2, depth+1, maxDepth);
    } else {
      Triangle sub;
      sub.x[0] = x0;
      sub.x[1] = x1;
      sub.x[2] = x2;
      sub.area = fabs(chi[0]*tau[1] - chi[1]*tau[0]);
      subTris.push_back(sub);
    }
  }
  
public:
  std::vector<Triangle> subTris;
  unsigned maxDepth;

  explicit TriangleRefiner() {}
  
  void refine(unsigned maxRefinementDepth) {
    subTris.clear();
    maxDepth = maxRefinementDepth;
    refine(glm::dvec2(0.0, 0.0), glm::dvec2(1.0, 0.0), glm::dvec2(0.0, 1.0), 0, maxRefinementDepth);
  }
};



#endif // TRIANGLE_REFINER_H
