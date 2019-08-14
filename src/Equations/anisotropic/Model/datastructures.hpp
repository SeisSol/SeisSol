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
        real c11;
        real c12;
        real c13;
        real c14;
        real c15;
        real c16;
        real c22;
        real c23;
        real c24;
        real c25;
        real c26;
        real c33;
        real c34;
        real c35;
        real c36;
        real c44;
        real c45;
        real c46;
        real c55;
        real c56;
        real c66;

        void getRotatedMaterialCoefficients( real      i_N[36],
                                             Material& o_Material ) 
        {
          o_Material.rho = rho;
          using Matrix66 = Eigen::Matrix<real, 6, 6>;
          Matrix66 N = Matrix66(i_N);
          Matrix66 C = Matrix66();
          C << c11, c12, c13, c14, c15, c16,
               c12, c22, c23, c24, c25, c26,
               c13, c23, c33, c34, c35, c36,
               c14, c24, c34, c44, c45, c46,
               c15, c25, c35, c45, c55, c56,
               c16, c26, c36, c46, c56, c66;
          Matrix66 rotatedC = N.transpose()*C*N;
          o_Material.c11 = rotatedC(0,0);
          o_Material.c12 = rotatedC(0,1);
          o_Material.c13 = rotatedC(0,2);
          o_Material.c14 = rotatedC(0,3);
          o_Material.c15 = rotatedC(0,4);
          o_Material.c16 = rotatedC(0,5);
          o_Material.c22 = rotatedC(1,1);
          o_Material.c23 = rotatedC(1,2);
          o_Material.c24 = rotatedC(1,3);
          o_Material.c25 = rotatedC(1,4);
          o_Material.c26 = rotatedC(1,5);
          o_Material.c33 = rotatedC(2,2);
          o_Material.c34 = rotatedC(2,3);
          o_Material.c35 = rotatedC(2,4);
          o_Material.c36 = rotatedC(2,5);
          o_Material.c44 = rotatedC(3,3);
          o_Material.c45 = rotatedC(3,4);
          o_Material.c46 = rotatedC(3,5);
          o_Material.c55 = rotatedC(4,4);
          o_Material.c56 = rotatedC(4,5);
          o_Material.c66 = rotatedC(5,5);
        }
    };
    struct LocalData {};
    struct NeighborData {};
  }
}

#endif
