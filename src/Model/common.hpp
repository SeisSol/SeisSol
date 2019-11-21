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
 
#ifndef MODEL_COMMON_HPP_
#define MODEL_COMMON_HPP_

#include <Eigen/Eigen>
#include <Initializer/typedefs.hpp>
#include <generated_code/init.h>
#include <Geometry/MeshTools.h>
#include <Numerical_aux/Transformation.h>
#include <Model/common_datastructures.hpp>
#include <Equations/elastic/Model/datastructures.hpp>
#include <Equations/anisotropic/Model/datastructures.hpp>
#include <Equations/viscoelastic2/Model/datastructures.hpp>
#include <iostream>

namespace seissol {
  namespace model {
    using Matrix99 = Eigen::Matrix<real, 9, 9>;

    template<typename Tmaterial, typename Tmatrix>
    void getTransposedCoefficientMatrix( Tmaterial const& i_material,
                                         unsigned         i_dim,
                                         Tmatrix&         o_M );
    
    template<typename Tmaterial, typename T>
    void getTransposedSourceCoefficientTensor(  Tmaterial const& material,
                                                T& E);

    template<typename Tmaterial, typename Tloc, typename Tneigh>
    void getTransposedGodunovState( Tmaterial const&  local,
                                    Tmaterial const&  neighbor,
                                    enum ::faceType   faceType,
                                    Tloc&             QgodLocal,
                                    Tneigh&           QgodNeighbor );

    template<typename T>
    void getTransposedFreeSurfaceGodunovState( T& QgodLocal,
                                               T& QgodNeighbor,
                                               Matrix99& R);

    template<typename T>
    void getPlaneWaveOperator( T const& material,
                               double const n[3],
                               std::complex<real> Mdata[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES] );

    template<typename T>
    void setMaterial( double* i_materialVal,
                      int i_numMaterialVals,
                      T* o_material ); 

    template<typename T>
    void initializeSpecificLocalData( T const&,
                                      LocalData* );

    template<typename T>
    void initializeSpecificNeighborData(  T const&,
                                          NeighborData* );

    void getBondMatrix( VrtxCoords const i_normal,
                        VrtxCoords const i_tangent1,
                        VrtxCoords const i_tangent2,
                        real* o_N );

    void getFaceRotationMatrix( VrtxCoords const i_normal,
                                VrtxCoords const i_tangent1,
                                VrtxCoords const i_tangent2,
                                init::T::view::type& o_T,
                                init::Tinv::view::type& o_Tinv );
  }
}

template<typename T>
void seissol::model::getPlaneWaveOperator(  T const& material,
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
}

template<typename T>
void seissol::model::getTransposedFreeSurfaceGodunovState( T&                         QgodLocal,
                                                           T&                         QgodNeighbor,
                                                           Eigen::Matrix<real, 9, 9>& R)
{
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      QgodNeighbor(i,j) = std::numeric_limits<double>::signaling_NaN();
    }
  }

  QgodLocal.setZero();
  std::array<int, 3> traction_indices = {0,3,5};
  std::array<int, 3> velocity_indices = {6,7,8};
  using Matrix33 = Eigen::Matrix<real, 3, 3>;
  Matrix33 R11 = R(traction_indices, {0,1,2});
  Matrix33 R21 = R(velocity_indices, {0,1,2});
  auto S = - (R21 * R11.inverse()).eval();

  //set lower left block
  int row = 0;
  for (auto &t: traction_indices) {
    int col = 0;
    for (auto &v: velocity_indices) {
      QgodLocal(t, v) = S(row, col);
      col++;
    }
    row++;
  }
  //set lower right block
  for (auto &v : velocity_indices) {
    QgodLocal(v, v) = 1.0;
  }
}


#endif
