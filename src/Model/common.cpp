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
 
#include <Model/common.hpp>
#include <cmath>

void seissol::model::getTransposedElasticCoefficientMatrix( seissol::model::ElasticMaterial const& i_material,
                                                            unsigned i_dim,
                                                            MatrixView o_M )
{
  o_M.setZero();

  real lambda2mu = i_material.lambda + 2.0 * i_material.mu;
  real rhoInv = 1.0 / i_material.rho;

  switch (i_dim)
  {
    case 0:
      o_M(6,0) = -lambda2mu;
      o_M(6,1) = -i_material.lambda;
      o_M(6,2) = -i_material.lambda;
      o_M(7,3) = -i_material.mu;
      o_M(8,5) = -i_material.mu;
      o_M(0,6) = -rhoInv;
      o_M(3,7) = -rhoInv;
      o_M(5,8) = -rhoInv;
      break;

    case 1:
      o_M(7,0) = -i_material.lambda;
      o_M(7,1) = -lambda2mu;
      o_M(7,2) = -i_material.lambda;
      o_M(6,3) = -i_material.mu;
      o_M(8,4) = -i_material.mu;
      o_M(3,6) = -rhoInv;
      o_M(1,7) = -rhoInv;
      o_M(4,8) = -rhoInv;
      break;

    case 2:
      o_M(8,0) = -i_material.lambda;
      o_M(8,1) = -i_material.lambda;
      o_M(8,2) = -lambda2mu;
      o_M(7,4) = -i_material.mu;
      o_M(6,5) = -i_material.mu;
      o_M(5,6) = -rhoInv;
      o_M(4,7) = -rhoInv;
      o_M(2,8) = -rhoInv;
      break;
      
    default:
      break;
  }
}

void seissol::model::getTransposedElasticGodunovState( seissol::model::ElasticMaterial const& local,
                                                       seissol::model::ElasticMaterial const& neighbor,
                                                       DenseMatrixView<9, 9> QgodLocal,
                                                       DenseMatrixView<9, 9> QgodNeighbor )
{
  QgodNeighbor.setZero();
  
  real cpL = sqrt((local.lambda + 2.0 * local.mu)       / local.rho);
  real cpN = sqrt((neighbor.lambda + 2.0 * neighbor.mu) / neighbor.rho);
  real csL = sqrt(local.mu / local.rho);
  real csN = sqrt(neighbor.mu / neighbor.rho);
  
  real constP = cpN * (local.lambda + 2.0 * local.mu) + cpL * (neighbor.lambda + 2.0 * neighbor.mu);
  real constS = csN * local.mu + csL * neighbor.mu;
  
  QgodNeighbor(0,0) = cpN * (local.lambda + 2.0 * local.mu) / constP;
  QgodNeighbor(6,0) = (local.lambda + 2.0 * local.mu) * (neighbor.lambda + 2.0 * neighbor.mu) / constP;
  QgodNeighbor(0,1) = cpN * local.lambda / constP;
  QgodNeighbor(6,1) = local.lambda * (neighbor.lambda + 2.0 * neighbor.mu) / constP;
  QgodNeighbor(0,2) = QgodNeighbor(0,1);
  QgodNeighbor(6,2) = QgodNeighbor(6,1);
  QgodNeighbor(3,3) = csN * local.mu / constS;
  QgodNeighbor(7,3) = local.mu * neighbor.mu / constS;
  QgodNeighbor(5,5) = QgodNeighbor(3,3);
  QgodNeighbor(8,5) = QgodNeighbor(7,3);
  QgodNeighbor(0,6) = cpL * cpN / constP;
  QgodNeighbor(6,6) = cpL * (neighbor.lambda + 2.0 * neighbor.mu) / constP;
  QgodNeighbor(3,7) = csL * csN / constS;
  QgodNeighbor(7,7) = csL * neighbor.mu / constS;
  QgodNeighbor(5,8) = QgodNeighbor(3,7);
  QgodNeighbor(8,8) = QgodNeighbor(7,7);

  // QgodLocal = I - QgodNeighbor
  for (unsigned idx = 0; idx < QgodLocal.rows() * QgodLocal.cols(); ++idx) {
    QgodLocal.data[idx] = -QgodNeighbor.data[idx];
  }  
  for (unsigned idx = 0; idx < QgodLocal.rows() && idx < QgodLocal.cols(); ++idx) {
    QgodLocal(idx, idx) += 1.0;
  }
}

void seissol::model::applyBoundaryConditionToElasticFluxSolver( enum ::faceType type,
                                                                DenseMatrixView<9, seissol::model::AminusT::cols> Fneighbor )
{
  if (type == freeSurface) {
    // Gamma is a diagonal matrix
    real Gamma[] = { -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0 };
    // Gamma^T * Fneighbor
    for (unsigned j = 0; j < Fneighbor.cols(); ++j) {
      for (unsigned i = 0; i < 9; ++i) {
        Fneighbor(i,j) *= Gamma[i];
      }
    }
  }
}
