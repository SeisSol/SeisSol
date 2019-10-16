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

        void getFullElasticTensor(real fullTensor[81]) {
          fullTensor[0]  = c(0,0);
          fullTensor[1]  = c(0,5);
          fullTensor[2]  = c(0,4);
          fullTensor[3]  = c(0,5);
          fullTensor[4]  = c(0,1);
          fullTensor[5]  = c(0,3);
          fullTensor[6]  = c(0,4);
          fullTensor[7]  = c(0,3);
          fullTensor[8]  = c(0,2);
          fullTensor[9]  = c(0,5);
          fullTensor[10] = c(5,5);
          fullTensor[11] = c(4,5);
          fullTensor[12] = c(5,5);
          fullTensor[13] = c(1,5);
          fullTensor[14] = c(3,5);
          fullTensor[15] = c(4,5);
          fullTensor[16] = c(3,5);
          fullTensor[17] = c(2,5);
          fullTensor[18] = c(0,4);
          fullTensor[19] = c(4,5);
          fullTensor[20] = c(4,4);
          fullTensor[21] = c(4,5);
          fullTensor[22] = c(1,4);
          fullTensor[23] = c(2,4);
          fullTensor[24] = c(4,4);
          fullTensor[25] = c(3,4);
          fullTensor[26] = c(2,4);
          fullTensor[27] = c(0,5);
          fullTensor[28] = c(5,5);
          fullTensor[29] = c(4,5);
          fullTensor[30] = c(5,5);
          fullTensor[31] = c(1,5);
          fullTensor[32] = c(3,5);
          fullTensor[33] = c(4,5);
          fullTensor[34] = c(3,5);
          fullTensor[35] = c(2,5);
          fullTensor[36] = c(0,1);
          fullTensor[37] = c(1,5);
          fullTensor[38] = c(1,4);
          fullTensor[39] = c(1,5);
          fullTensor[40] = c(1,1);
          fullTensor[41] = c(1,3);
          fullTensor[42] = c(1,4);
          fullTensor[43] = c(1,3);
          fullTensor[44] = c(1,2);
          fullTensor[45] = c(0,3);
          fullTensor[46] = c(3,5);
          fullTensor[47] = c(3,4);
          fullTensor[48] = c(3,5);
          fullTensor[49] = c(1,3);
          fullTensor[50] = c(3,3);
          fullTensor[51] = c(3,4);
          fullTensor[52] = c(3,3);
          fullTensor[53] = c(2,3);
          fullTensor[54] = c(0,4);
          fullTensor[55] = c(4,5);
          fullTensor[56] = c(4,4);
          fullTensor[57] = c(4,5);
          fullTensor[58] = c(1,4);
          fullTensor[59] = c(3,4);
          fullTensor[60] = c(4,4);
          fullTensor[61] = c(3,4);
          fullTensor[62] = c(2,4);
          fullTensor[63] = c(0,3);
          fullTensor[64] = c(3,5);
          fullTensor[65] = c(3,4);
          fullTensor[66] = c(3,5);
          fullTensor[67] = c(1,3);
          fullTensor[68] = c(3,3);
          fullTensor[69] = c(3,4);
          fullTensor[70] = c(3,3);
          fullTensor[71] = c(2,3);
          fullTensor[72] = c(0,3);
          fullTensor[73] = c(2,5);
          fullTensor[74] = c(2,4);
          fullTensor[75] = c(2,5);
          fullTensor[76] = c(1,2);
          fullTensor[77] = c(2,3);
          fullTensor[78] = c(2,4);
          fullTensor[79] = c(2,3);
          fullTensor[80] = c(2,2);
        }
    };
    struct LocalData {};
    struct NeighborData {};
  }
}

#endif
