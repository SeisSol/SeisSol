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

#include <cmath>

#include <Model/Setup.h>

#include <Model/common.hpp>
#include <Kernels/common.hpp>
#include <Numerical_aux/Transformation.h>

#include <yateto/TensorView.h>
#include <generated_code/init.h>
#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <stdexcept>
#include <iostream>

template<typename T>
void getTransposedCoefficientMatrix( seissol::model::Material const& material,
                                     unsigned dim,
                                     T& AT)
{
  switch(dim){
    case 0:
      AT(0,6) = -1/(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1));
      AT(0,10) = -1/(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity);
      AT(3,7) = -1/(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1));
      AT(3,11) = -1/(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity);
      AT(5,8) = -1/(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1));
      AT(5,12) = -1/(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity);
      AT(6,0) = -1.0/3.0*material.bulk_solid*pow(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid, 2)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid) - material.lambda - 2*material.mu;
      AT(6,1) = -1.0/3.0*material.bulk_solid*pow(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid, 2)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid) - material.lambda;
      AT(6,2) = -1.0/3.0*material.bulk_solid*pow(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid, 2)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid) - material.lambda;
      AT(6,9) = -material.bulk_solid*(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid);
      AT(7,3) = -material.mu;
      AT(8,5) = -material.mu;
      AT(9,6) = -material.porosity/(material.tortuosity*(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1)));
      AT(9,10) = -(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/(material.rho_fluid*(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity));
      AT(10,0) = material.bulk_solid*(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid);
      AT(10,1) = material.bulk_solid*(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid);
      AT(10,2) = material.bulk_solid*(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid);
      AT(10,9) = 3*material.bulk_solid/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid);
      break;
    case 1:
      AT(1,7) = -1/(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1));
      AT(1,11) = -1/(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity);
      AT(3,6) = -1/(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1));
      AT(3,10) = -1/(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity);
      AT(4,8) = -1/(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1));
      AT(4,12) = -1/(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity);
      AT(6,3) = -material.mu;
      AT(7,0) = -1.0/3.0*material.bulk_solid*pow(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid, 2)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid) - material.lambda;
      AT(7,1) = -1.0/3.0*material.bulk_solid*pow(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid, 2)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid) - material.lambda - 2*material.mu;
      AT(7,2) = -1.0/3.0*material.bulk_solid*pow(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid, 2)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid) - material.lambda;
      AT(7,9) = -material.bulk_solid*(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid);
      AT(8,4) = -material.mu;
      AT(9,7) = -material.porosity/(material.tortuosity*(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1)));
      AT(9,11) = -(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/(material.rho_fluid*(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity));
      AT(11,0) = material.bulk_solid*(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid);
      AT(11,1) = material.bulk_solid*(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid);
      AT(11,2) = material.bulk_solid*(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid);
      AT(11,9) = 3*material.bulk_solid/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid);
      break;
    case 2:
      AT(2,8) = -1/(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1));
      AT(2,12) = -1/(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity);
      AT(4,7) = -1/(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1));
      AT(4,11) = -1/(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity);
      AT(5,6) = -1/(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1));
      AT(5,10) = -1/(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity);
      AT(6,5) = -material.mu;
      AT(7,4) = -material.mu;
      AT(8,0) = -1.0/3.0*material.bulk_solid*pow(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid, 2)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid) - material.lambda;
      AT(8,1) = -1.0/3.0*material.bulk_solid*pow(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid, 2)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid) - material.lambda;
      AT(8,2) = -1.0/3.0*material.bulk_solid*pow(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid, 2)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid) - material.lambda - 2*material.mu;
      AT(8,9) = -material.bulk_solid*(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid);
      AT(9,8) = -material.porosity/(material.tortuosity*(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1)));
      AT(9,12) = -(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/(material.rho_fluid*(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity));
      AT(12,0) = material.bulk_solid*(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid);
      AT(12,1) = material.bulk_solid*(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid);
      AT(12,2) = material.bulk_solid*(-3 + (3*material.lambda + 2*material.mu)/material.bulk_solid)/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid);
      AT(12,9) = 3*material.bulk_solid/(3*material.porosity*(-1 + material.bulk_solid/material.bulk_fluid) + 3 - (3*material.lambda + 2*material.mu)/material.bulk_solid);
      break;
    default:
      break;
  }

}

void seissol::model::getTransposedCoefficientMatrix( Material const& material,
                                                     unsigned        dim,
                                                     init::star::view<0>::type& AT )
{
  AT.setZero();
  ::getTransposedCoefficientMatrix(material, dim, AT);
}

template<typename T>
void getTransposedSourceCoefficientTensor(seissol::model::Material const& material, T& ET) {
  ET.setZero();
  ET(10,10) = material.viscosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/(material.permeability*material.rho_fluid*(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity));
  ET(11,11) = material.viscosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/(material.permeability*material.rho_fluid*(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity));
  ET(12,12) = material.viscosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/(material.permeability*material.rho_fluid*(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity));

  ET(10,6) = material.porosity*material.viscosity/(material.permeability*material.tortuosity*(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1)));
  ET(11,7) = material.porosity*material.viscosity/(material.permeability*material.tortuosity*(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1)));
  ET(12,8) = material.porosity*material.viscosity/(material.permeability*material.tortuosity*(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1)));
}

void seissol::model::getPlaneWaveOperator(  Material const& material,
                                            double const n[3],
                                            std::complex<real> Mdata[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES] )
{
  yateto::DenseTensorView<2,std::complex<real>> M(Mdata, {NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES});
  M.setZero();

  real data[NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES];
  yateto::DenseTensorView<2,real> Coeff(data, {NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES});

  for (unsigned d = 0; d < 3; ++d) {
    Coeff.setZero();
    ::getTransposedCoefficientMatrix(material, d, Coeff);

    for (unsigned i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
      for (unsigned j = 0; j < NUMBER_OF_QUANTITIES; ++j) {
        M(i,j) += n[d] * Coeff(j,i);
      }
    }
  }

  ::getTransposedSourceCoefficientTensor(material, Coeff);

  for (unsigned i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
    for (unsigned j = 0; j < NUMBER_OF_QUANTITIES; ++j) {
      M(i,j) -= std::complex<real>(0.0, Coeff(j,i));
    }
  }
}

void seissol::model::getTransposedGodunovState( Material const&                   local,
                                                Material const&                   neighbor,
                                                FaceType                          faceType,
                                                init::QgodLocal::view::type&      QgodLocal,
                                                init::QgodNeighbor::view::type&   QgodNeighbor )
{
  constexpr auto tolerance = 1.0e-4;

  using Matrix = Eigen::Matrix<double, 13, 13, Eigen::ColMajor>;
  using CMatrix = Eigen::Matrix<std::complex<double>, 13, 13, Eigen::ColMajor>;
  auto eigenDecomposition = [&tolerance](Material const& material) {
    Matrix t = Matrix::Zero();
    ::getTransposedCoefficientMatrix(material, 0, t);
    Eigen::EigenSolver<Matrix> es;
    es.compute(t.transpose());

#ifndef NDEBUG
    auto evs = es.eigenvalues();
    int ev_neg = 0;
    int ev_pos = 0;
    for (int i = 0; i < 13; ++i) {
      if (evs(i).real() < -tolerance) {
        ++ev_neg;
      } else if (evs(i).real() > tolerance) {
        ++ev_pos;
      }
    }
    assert(ev_neg == 4);
    assert(ev_pos == 4);
#endif
    return es;
  };
  auto eigen_local = eigenDecomposition(local);
  auto eigen_neighbor = eigenDecomposition(neighbor); 
  Matrix chi_minus = Matrix::Zero();
  Matrix chi_plus = Matrix::Zero();
  for(int i = 0; i < 13; i++) {
    if(eigen_local.eigenvalues()[i].real() < tolerance) {
      chi_minus(i,i) = 1.0;
    }
    if(eigen_neighbor.eigenvalues()[i].real() > tolerance) {
      chi_plus(i,i) = 1.0;
    }
  }
  CMatrix R;
  R = eigen_local.eigenvectors() * chi_minus + eigen_neighbor.eigenvectors() * chi_plus;
  //R = eigen_local.eigenvectors();
  if (faceType == FaceType::freeSurface) {
    logWarning() << "Poroelastic Free Surface is not tested yet.";
    Matrix R_real = R.real().eval();
    getTransposedFreeSurfaceGodunovState(local, QgodLocal, QgodNeighbor, R_real);
  } else {
    CMatrix godunov_minus = ((R*chi_minus)*R.inverse()).eval();
    CMatrix godunov_plus = ((R*chi_plus)*R.inverse()).eval();


    for (unsigned i = 0; i < QgodLocal.shape(1); ++i) {
      for (unsigned j = 0; j < QgodLocal.shape(0); ++j) {
        QgodLocal(i,j) = godunov_plus(j,i).real();
        QgodNeighbor(i,j) = godunov_minus(j,i).real();
      }
    }
  }

  /*QgodLocal.setZero();
  QgodNeighbor.setZero();
  for (unsigned i = 0; i < 13; ++i) {
    QgodLocal(i,i) = 0.5;
    QgodNeighbor(i,i) = 0.5;
  }*/

}

void seissol::model::setMaterial( double* i_materialVal,
                                  int i_numMaterialVals,
                                  seissol::model::Material* o_material )
{
  assert(i_numMaterialVals == 10);
  o_material->bulk_solid = i_materialVal[0];
  o_material->rho = i_materialVal[1]; 
  o_material->lambda = i_materialVal[2];    
  o_material->mu = i_materialVal[3];
  o_material->porosity = i_materialVal[4]; 
  o_material->permeability = i_materialVal[5];
  o_material->tortuosity = i_materialVal[6];
  o_material->bulk_fluid = i_materialVal[7];
  o_material->rho_fluid = i_materialVal[8];
  o_material->viscosity = i_materialVal[9];  
}


void seissol::model::getFaceRotationMatrix( VrtxCoords const i_normal,
                                            VrtxCoords const i_tangent1,
                                            VrtxCoords const i_tangent2,
                                            init::T::view::type& o_T,
                                            init::Tinv::view::type& o_Tinv )
{
  o_T.setZero();
  o_Tinv.setZero();
  
  seissol::transformations::symmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, o_T);
  seissol::transformations::inverseSymmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, o_Tinv);
  
  seissol::transformations::tensor1RotationMatrix(i_normal, i_tangent1, i_tangent2, o_T, 6, 6);
  seissol::transformations::inverseTensor1RotationMatrix(i_normal, i_tangent1, i_tangent2, o_Tinv, 6, 6);
  
  unsigned origin = 10; 
  seissol::transformations::tensor1RotationMatrix(i_normal, i_tangent1, i_tangent2, o_T, origin, origin);
  seissol::transformations::inverseTensor1RotationMatrix(i_normal, i_tangent1, i_tangent2, o_Tinv, origin, origin);

  o_T(9, 9) = 1;
  o_Tinv(9,9) = 1;
}

void seissol::model::initializeSpecificLocalData( seissol::model::Material const& material,
                                                  seissol::model::LocalData* localData )
{
  auto sourceMatrix = init::ET::view::create(localData->sourceMatrix);
  sourceMatrix.setZero();
  getTransposedSourceCoefficientTensor(material, sourceMatrix);

  using Matrix = Eigen::Matrix<real, CONVERGENCE_ORDER, CONVERGENCE_ORDER>;
  using Vector = Eigen::Matrix<real, CONVERGENCE_ORDER, 1>;

  auto Zinv = init::Zinv::view::create(localData->Zinv); 
  for(int i = 0; i < NUMBER_OF_QUANTITIES; i++) {
    Matrix Z(init::Z::Values);
    if(i >= 10) {
      for(int j = 0; j < CONVERGENCE_ORDER; j++) {
        Z(j,j) = Z(j,j) - sourceMatrix(i,i);
      }
    }
    auto solver = Z.colPivHouseholderQr();
    for(int col = 0; col < CONVERGENCE_ORDER; col++) {
      Vector rhs = Vector::Zero();
      rhs(col) = 1.0;
      auto Zinv_col = solver.solve(rhs);
        for(int row = 0; row < CONVERGENCE_ORDER; row++) {
          Zinv(i,row,col) = Zinv_col(row);
      }
    }
  }

  
}

void seissol::model::initializeSpecificNeighborData(  seissol::model::Material const&,
                                                      seissol::model::Material const (&)[4],
                                                      seissol::model::NeighborData* )
{
}
