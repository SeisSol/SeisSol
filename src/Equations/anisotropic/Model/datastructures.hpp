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

#ifndef MODEL_DATASTRUCTURES_H_
#define MODEL_DATASTRUCTURES_H_

#include <Model/common_datastructures.hpp>
#include <Eigen/Eigen>

namespace seissol {
  namespace model {
    struct Material {
        real rho;
        real c_store[21];
        int c_addr(unsigned i, unsigned j) const {
          auto row = std::min(i,j); 
          auto col = std::max(i,j); 
          auto offset = row*(11-row)/2;
          return offset + col;
        }
        real c(unsigned i, unsigned j) const {
          return c_store[c_addr(i,j)];
        }

        void getRotatedMaterialCoefficients( real      i_N[36],
                                             Material& o_Material ) 
        {
          o_Material.rho = rho;
          using Matrix66 = Eigen::Matrix<real, 6, 6>;
          Matrix66 N = Matrix66(i_N);
          Matrix66 C = Matrix66();
          for(int i = 0; i < 6; i++)
            for(int j = 0; j < 6; j++)
              C(i,j) = c(i,j);
          Matrix66 rotatedC = N.transpose()*C*N;
          for(int i = 0; i < 6; i++)
            for(int j = i; j < 6; j++)
              o_Material.c_store[c_addr(i,j)] = rotatedC(i,j);
        }
    };
    struct LocalData {};
    struct NeighborData {};
  }
}

#endif
