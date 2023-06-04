/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT in.tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015 - 2019, SeisSol Group
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

#ifndef VISCOELASTIC_SETUP_H_
#define VISCOELASTIC_SETUP_H_

#include <Model/common.hpp>
#include <Kernels/common.hpp>
#include <Numerical_aux/Transformation.h>

#include <yateto/TensorView.h>
#include <generated_code/init.h>

namespace seissol {
  namespace model {
    template<typename T>
    inline void getTransposedViscoelasticCoefficientMatrix( real            i_omega,
                                                            unsigned        i_dim,
                                                            unsigned        mech,
                                                            T&              o_M )
    {
      unsigned col = 9 + mech * 6;
      switch (i_dim)
      {
        case 0:
          o_M(6, col)     = -i_omega;
          o_M(7, col + 3) = -0.5 * i_omega;
          o_M(8, col + 5) = -0.5 * i_omega;
          break;

        case 1:
          o_M(7, col + 1) = -i_omega;
          o_M(6, col + 3) = -0.5 * i_omega;
          o_M(8, col + 4) = -0.5 * i_omega;
          break;

        case 2:
          o_M(8, col + 2) = -i_omega;
          o_M(7, col + 4) = -0.5 * i_omega;
          o_M(6, col + 5) = -0.5 * i_omega;
          break;
      }
    }

    template<typename T, std::size_t Mechanisms>
    inline void getTransposedSourceCoefficientTensor( ViscoElasticMaterial<Mechanisms> const& material, T& sourceMatrix )
    {
      sourceMatrix.setZero();
    
      //       | E_1^T |
      // E^T = |  ...  |
      //       | E_L^T |
      for (unsigned mech = 0; mech < Mechanisms; ++mech) {
        unsigned offset = ViscoElasticMaterial<Mechanisms>::NumberOfQuantities + mech * ViscoElasticMaterial<Mechanisms>::NumberPerMechanism;
        double const* theta = material.theta[mech];
        sourceMatrix(offset,     0) = theta[0];
        sourceMatrix(offset + 1, 0) = theta[1];
        sourceMatrix(offset + 2, 0) = theta[1];
        sourceMatrix(offset,     1) = theta[1];
        sourceMatrix(offset + 1, 1) = theta[0];
        sourceMatrix(offset + 2, 1) = theta[1];  
        sourceMatrix(offset,     2) = theta[1];
        sourceMatrix(offset + 1, 2) = theta[1];
        sourceMatrix(offset + 2, 2) = theta[0];  
        sourceMatrix(offset + 3, 3) = theta[2];
        sourceMatrix(offset + 4, 4) = theta[2];
        sourceMatrix(offset + 5, 5) = theta[2];    
      }
      
      // E' = diag(-omega_1 I, ..., -omega_L I)
      for (unsigned mech = 0; mech < Mechanisms; ++mech) {
        for (unsigned i = 0; i < ViscoElasticMaterial<Mechanisms>::NumberPerMechanism; ++i) {
          unsigned idx = ViscoElasticMaterial<Mechanisms>::NumberOfQuantities + mech * ViscoElasticMaterial<Mechanisms>::NumberPerMechanism + i;
          sourceMatrix(idx, idx) = -material.omega[mech];
        }
      }
    }
    
    template<typename T, std::size_t Mechanisms>
    inline void getTransposedCoefficientMatrix(ViscoElasticMaterial<Mechanisms> const& i_material, unsigned i_dim, T& AT)
    {
      getTransposedCoefficientMatrix(dynamic_cast<ElasticMaterial const&>(i_material), i_dim, AT);
    
      for (unsigned mech = 0; mech < Mechanisms; ++mech) {
        getTransposedViscoelasticCoefficientMatrix( i_material.omega[mech],
                                                    i_dim,
                                                    mech,
                                                    AT );
      }
    }
    
    template<std::size_t Mechanisms>
    inline void getTransposedGodunovState( ViscoElasticMaterial<Mechanisms> const& local,
                                           ViscoElasticMaterial<Mechanisms> const& neighbor,
                                           FaceType faceType,
                                           init::QgodLocal::view::type& QgodLocal,
                                           init::QgodNeighbor::view::type& QgodNeighbor)
    {
      getTransposedGodunovState( dynamic_cast<ElasticMaterial const&>(local),
                                 dynamic_cast<ElasticMaterial const&>(neighbor), 
                                 faceType, 
                                 QgodLocal, 
                                 QgodNeighbor);
    }

    template<std::size_t Mechanisms>
    inline void initializeSpecificLocalData( ViscoElasticMaterial<Mechanisms> const& material,
                                             real timeStepWidth,
                                             ViscoElasticLocalData* localData )
    {
      auto sourceMatrix = init::ET::view::create(localData->sourceMatrix);
      getTransposedSourceCoefficientTensor(material, sourceMatrix);
    }
  }
}

#endif
