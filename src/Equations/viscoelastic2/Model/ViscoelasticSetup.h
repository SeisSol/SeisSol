/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT in.tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015 - 2020, SeisSol Group
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

#ifndef VISCOELASTIC2_SETUP_H_
#define VISCOELASTIC2_SETUP_H_


#include "Model/common.hpp"
#include "Kernels/common.hpp"
#include "Numerical_aux/Transformation.h"

#include <TensorForge/TensorView.h>
#include "generated_code/init.h"

namespace seissol {
  namespace model {
    using Matrix99 = Eigen::Matrix<double, 9, 9>;

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

/*
 * The new implemenation of attenuation (viscoelastic2) is considered standard. This part will be 
 * used unless the old attenuation (viscoelastic) implementation is chosen. 
 */
    template<typename T>
    void getTransposedSourceCoefficientTensor(  ViscoElasticMaterial const& material,
                                                T& E )
      {
        for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
          double const* theta = material.theta[mech];
          E(0, mech, 0) = theta[0];
          E(1, mech, 0) = theta[1];
          E(2, mech, 0) = theta[1];
          E(0, mech, 1) = theta[1];
          E(1, mech, 1) = theta[0];
          E(2, mech, 1) = theta[1];
          E(0, mech, 2) = theta[1];
          E(1, mech, 2) = theta[1];
          E(2, mech, 2) = theta[0];
          E(3, mech, 3) = theta[2];
          E(4, mech, 4) = theta[2];
          E(5, mech, 5) = theta[2];
        }
      }

    template<typename T>
    void getTransposedCoefficientMatrix( ViscoElasticMaterial const&    i_material,
                                         unsigned                       i_dim,
                                         T&                             AT )
      {
        getTransposedCoefficientMatrix(dynamic_cast<ElasticMaterial const&>(i_material), i_dim, AT);

        getTransposedViscoelasticCoefficientMatrix( 1.0, i_dim, 0, AT );
      }


    template<>
    inline void getTransposedGodunovState( ViscoElasticMaterial const&       local,
                                    ViscoElasticMaterial const&       neighbor,
                                    FaceType                          faceType,
                                    init::QgodLocal::view::type&      QgodLocal,
                                    init::QgodNeighbor::view::type&   QgodNeighbor )
    {
      getTransposedGodunovState( dynamic_cast<ElasticMaterial const&>(local),
                                 dynamic_cast<ElasticMaterial const&>(neighbor), 
                                 faceType, 
                                 QgodLocal, 
                                 QgodNeighbor);
    }

    template<>
    inline void getPlaneWaveOperator( ViscoElasticMaterial const& material,
                                      double const n[3],
                                      std::complex<double> Mdata[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES] )
      {
        yateto::DenseTensorView<2,std::complex<double>> M(Mdata, {NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES});
        M.setZero();

        double data[NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES];
        yateto::DenseTensorView<2,double> Coeff(data, {NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES});

        for (unsigned d = 0; d < 3; ++d) {
          Coeff.setZero();
          getTransposedCoefficientMatrix(material, d, Coeff);
          for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
            getTransposedViscoelasticCoefficientMatrix( material.omega[mech],
                d,
                mech,
                Coeff );
          }

          for (unsigned i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
            for (unsigned j = 0; j < NUMBER_OF_QUANTITIES; ++j) {
              M(i,j) += n[d] * Coeff(j,i);
            }
          }
        }
        double Edata[NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES];
        yateto::DenseTensorView<3,double> E(Edata, tensor::E::Shape);
        E.setZero();
        getTransposedSourceCoefficientTensor(material, E);
        Coeff.setZero();
        for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
          unsigned offset = 9 + mech * 6;
          for (unsigned i = 0; i < tensor::E::Shape[0]; ++i) {
            for (unsigned j = 0; j < tensor::E::Shape[2]; ++j) {
              Coeff(offset + i, j) = E(i, mech, j);
            }
          }
        }

        // E' = diag(-omega_1 I, ..., -omega_L I)
        for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
          unsigned offset = 9 + 6*mech;
          yateto::DenseTensorView<2,double> ETblock(data + offset + offset * NUMBER_OF_QUANTITIES, {NUMBER_OF_QUANTITIES, 6});
          for (unsigned i = 0; i < 6; ++i) {
            ETblock(i, i) = -material.omega[mech];
          }
        }

        for (unsigned i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
          for (unsigned j = 0; j < NUMBER_OF_QUANTITIES; ++j) {
            M(i,j) -= std::complex<double>(0.0, Coeff(j,i));
          }
        }
      } 

    template<>
    inline void initializeSpecificLocalData( ViscoElasticMaterial const& material,
                                             real timeStepWidth,
                                             ViscoElasticLocalData* localData )
    {
      auto E = init::E::view::create(localData->E);
      E.setZero();
      getTransposedSourceCoefficientTensor(material, E);

      auto w = init::w::view::create(localData->w);
      auto W = init::W::view::create(localData->W);
      W.setZero();
      for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
        w(mech) = material.omega[mech];
        W(mech,mech) = -material.omega[mech];
      }
    }

    template<>
    inline void initializeSpecificNeighborData(  ViscoElasticMaterial const& localMaterial,
                                                 ViscoElasticNeighborData* neighborData )
    {
      // We only need the local omegas
      auto w = init::w::view::create(neighborData->w);
      for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
        w(mech) = localMaterial.omega[mech];
      }
    }
  }
}

#endif
