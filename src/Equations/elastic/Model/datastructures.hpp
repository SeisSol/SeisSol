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

#ifndef MODEL_ELASTIC_DATASTRUCTURES_H_
#define MODEL_ELASTIC_DATASTRUCTURES_H_

#include <Model/common_datastructures.hpp>
#include <cmath>

namespace seissol {
  namespace model {
    struct ElasticMaterial : Material {
      real lambda;
      real mu;

      virtual ~ElasticMaterial() {};

      Material getRotatedMaterialCoefficients(real i_N[36]) {
        // isotropic materials are rotationally invariant
        return *this;
      }

      void getFullElasticTensor(real fullTensor[81]) {
        fullTensor[0]  = lambda + 2*mu;
        fullTensor[1]  = 0.0;
        fullTensor[2]  = 0.0;
        fullTensor[3]  = 0.0;
        fullTensor[4]  = lambda;
        fullTensor[5]  = 0.0;
        fullTensor[6]  = 0.0;
        fullTensor[7]  = 0.0;
        fullTensor[8]  = lambda;
        fullTensor[9]  = 0.0;
        fullTensor[10] = 0.0;
        fullTensor[11] = 0.0;
        fullTensor[12] = 0.0;
        fullTensor[13] = 0.0;
        fullTensor[14] = 0.0;
        fullTensor[15] = 0.0;
        fullTensor[16] = 0.0;
        fullTensor[17] = 0.0;
        fullTensor[18] = 0.0;
        fullTensor[19] = 0.0;
        fullTensor[20] = mu;
        fullTensor[21] = 0.0;
        fullTensor[22] = 0.0;
        fullTensor[23] = 0.0;
        fullTensor[24] = mu;
        fullTensor[25] = 0.0;
        fullTensor[26] = 0.0;
        fullTensor[27] = 0.0;
        fullTensor[28] = 0.0;
        fullTensor[29] = 0.0;
        fullTensor[30] = 0.0;
        fullTensor[31] = 0.0;
        fullTensor[32] = 0.0;
        fullTensor[33] = 0.0;
        fullTensor[34] = 0.0;
        fullTensor[35] = 0.0;
        fullTensor[36] = lambda;
        fullTensor[37] = 0.0;
        fullTensor[38] = 0.0;
        fullTensor[39] = 0.0;
        fullTensor[40] = lambda + 2*mu;
        fullTensor[41] = 0.0;
        fullTensor[42] = 0.0;
        fullTensor[43] = 0.0;
        fullTensor[44] = lambda;
        fullTensor[45] = 0.0;
        fullTensor[46] = 0.0;
        fullTensor[47] = 0.0;
        fullTensor[48] = 0.0;
        fullTensor[49] = 0.0;
        fullTensor[50] = mu;
        fullTensor[51] = 0.0;
        fullTensor[52] = mu;
        fullTensor[53] = 0.0;
        fullTensor[54] = 0.0;
        fullTensor[55] = 0.0;
        fullTensor[56] = mu;
        fullTensor[57] = 0.0;
        fullTensor[58] = 0.0;
        fullTensor[59] = 0.0;
        fullTensor[60] = mu;
        fullTensor[61] = 0.0;
        fullTensor[62] = 0.0;
        fullTensor[63] = 0.0;
        fullTensor[64] = 0.0;
        fullTensor[65] = 0.0;
        fullTensor[66] = 0.0;
        fullTensor[67] = 0.0;
        fullTensor[68] = mu;
        fullTensor[69] = 0.0;
        fullTensor[70] = mu;
        fullTensor[71] = 0.0;
        fullTensor[72] = 0.0;
        fullTensor[73] = 0.0;
        fullTensor[74] = 0.0;
        fullTensor[75] = 0.0;
        fullTensor[76] = lambda;
        fullTensor[77] = 0.0;
        fullTensor[78] = 0.0;
        fullTensor[79] = 0.0;
        fullTensor[80] = lambda + 2*mu;
      }

      real getMaxWaveSpeed() {
        return getPWaveSpeed();
      }

      real getPWaveSpeed() {
        return std::sqrt((lambda + 2*mu) / rho);
      }

      real getSWaveSpeed() {
        return std::sqrt(mu / rho);
      }

      materialType getMaterialType() const {
        return elastic;
      }
    };

    struct ElastoPlasticMaterial : ElasticMaterial, Plasticity {
      virtual ~ElastoPlasticMaterial() {};
      materialType getMaterialType() const {
        return elastoplastic;
      }
    };
  }
}

#endif
