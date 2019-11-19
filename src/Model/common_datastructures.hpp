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
 
#ifndef MODEL_COMMONDATASTRUCTURES_HPP_
#define MODEL_COMMONDATASTRUCTURES_HPP_

#include <Kernels/precision.hpp>


namespace seissol {
  namespace model {
    enum materialType {
      elastic,
      viscoelastic,
      elastoplastic,
      viscoplastic,
      anisotropic
    };

    struct Material {
      real rho;
      virtual real getMaxWaveSpeed() {};
      virtual real getPWaveSpeed() {};
      virtual real getSWaveSpeed() {};
      virtual Material getRotatedMaterialCoefficients(real i_N[36]) {} ;
      virtual void getFullElasticTensor(real fullTensor[81]) {}; 
      virtual materialType getMaterialType() const {};
      virtual ~Material() {};
    };

    struct Plasticity {
      real bulkFriction;
      real plastCo;  
      real s_xx;   
      real s_yy;   
      real s_zz;   
      real s_xy;      
      real s_yz;      
      real s_xz;        
    };

    struct IsotropicWaveSpeeds {
      real density;
      real pWaveVelocity;
      real sWaveVelocity;
    };
#ifndef USE_VISCOELASTIC
    struct LocalData {
    };
    struct NeighborData {
    };
#endif
  }
  
}

#endif
