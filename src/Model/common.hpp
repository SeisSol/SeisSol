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

#include "utils/logger.h"
#include "Initializer/typedefs.hpp"
#include "generated_code/init.h"

namespace seissol {
  namespace model {
    using Matrix99 = Eigen::Matrix<double, 9, 9>;

    bool testIfAcoustic(real mu);

    template<typename T>
    void getTransposedElasticCoefficientMatrix(ElasticMaterial const& i_material,
                                               unsigned i_dim,
                                               T& o_M);

    template<typename Tloc, typename Tneigh>
    void getTransposedElasticGodunovState(Material const& local,
                                          Material const& neighbor,
                                          ::FaceType faceType,
                                          Tloc& QgodLocal,
                                          Tneigh& QgodNeighbor);

    template<typename T>
    void getTransposedFreeSurfaceGodunovState(const Material& local,
                                              T& QgodLocal,
                                              T& QgodNeighbor,
                                              Matrix99& R);
    template<typename T>
    void applyBoundaryConditionToElasticFluxSolver(::FaceType type,
                                                   T& QgodNeighbor);
  }
}

template<typename T>
void seissol::model::getTransposedElasticCoefficientMatrix(seissol::model::ElasticMaterial const& i_material,
                                                           unsigned i_dim,
                                                           T& o_M) {
  o_M.setZero();

  const real lambda2mu = i_material.lambda + 2.0 * i_material.mu;
  const real rhoInv = 1.0 / i_material.rho;

  switch (i_dim)
  {
    case 0:
      o_M(6,0) = -lambda2mu;
      o_M(6,1) = -i_material.lambda;
      o_M(6,2) = -i_material.lambda;
      o_M(7,3) = -i_material.mu;
      o_M(8,5) = -i_material.mu;
      o_M(0,6) = -rhoInv;
      if (!testIfAcoustic(i_material.mu)) {
        o_M(3,7) = -rhoInv;
        o_M(5,8) = -rhoInv;
      }
      break;

    case 1:
      o_M(7,0) = -i_material.lambda;
      o_M(7,1) = -lambda2mu;
      o_M(7,2) = -i_material.lambda;
      o_M(6,3) = -i_material.mu;
      o_M(8,4) = -i_material.mu;
      o_M(1,7) = -rhoInv;
      if (!testIfAcoustic(i_material.mu)) {
        o_M(3,6) = -rhoInv;
        o_M(4,8) = -rhoInv;
      }
      break;

    case 2:
      o_M(8,0) = -i_material.lambda;
      o_M(8,1) = -i_material.lambda;
      o_M(8,2) = -lambda2mu;
      o_M(7,4) = -i_material.mu;
      o_M(6,5) = -i_material.mu;
      o_M(2,8) = -rhoInv;
      if (!testIfAcoustic(i_material.mu)) {
        o_M(5,6) = -rhoInv;
        o_M(4,7) = -rhoInv;
      }
      break;
      
    default:
      break;
  }
}

template<typename T>
void seissol::model::getTransposedFreeSurfaceGodunovState(const Material& material,
                                                          T& QgodLocal,
                                                          T& QgodNeighbor,
                                                          Matrix99& R) {
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      QgodNeighbor(i,j) = std::numeric_limits<double>::signaling_NaN();
    }
  }

  QgodLocal.setZero();
  if (testIfAcoustic(material.mu)) {
    // Acoustic material only has one traction (=pressure) and one velocity comp.
    // relevant to the Riemann problem
    QgodLocal(0, 6) = -1 * R(6,0) * 1/R(0,0); // S
    QgodLocal(6, 6) = 1.0;
  } else {
    std::array<int, 3> traction_indices = {0,3,5};
    std::array<int, 3> velocity_indices = {6,7,8};
    using Matrix33 = Eigen::Matrix<double, 3, 3>;
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
}


template<typename Tloc, typename Tneigh>
void seissol::model::getTransposedElasticGodunovState(Material const& local,
                                                      Material const& neighbor,
                                                      FaceType faceType,
                                                      Tloc& QgodLocal,
                                                      Tneigh& QgodNeighbor) {
  QgodNeighbor.setZero();

  // Eigenvectors are precomputed
  Matrix99 R = Matrix99::Zero();

  if (testIfAcoustic(local.mu)) {
    R(0,0) = local.lambda;
    R(1,0) = local.lambda;
    R(2,0) = local.lambda;
    R(6,0) = std::sqrt((local.lambda) / local.rho);

    R(3,1) = 1.0;
    R(5,2) = 1.0;
  } else {
    R(0,0) = local.lambda + 2*local.mu;
    R(1,0) = local.lambda;
    R(2,0) = local.lambda;
    R(6,0) = std::sqrt((local.lambda + 2 * local.mu) / local.rho);

    R(3,1) = local.mu;
    R(7,1) = std::sqrt(local.mu / local.rho);

    R(5,2) = local.mu;
    R(8,2) = std::sqrt(local.mu / local.rho);
  }

  R(4,3) = 1;

  R(1,4) = 1;

  R(2,5) = 1;

  if (testIfAcoustic(neighbor.mu)) {
    R(7,6) = 1.0;
    R(8, 7) = 1.0;

    R(0,8) = neighbor.lambda;
    R(1,8) = neighbor.lambda;
    R(2,8) = neighbor.lambda;
    R(6,8) = -std::sqrt((neighbor.lambda + 2 * neighbor.mu) / neighbor.rho);
  } else {
    R(5,6) = neighbor.mu;
    R(8,6) = -std::sqrt(neighbor.mu / neighbor.rho);

    R(3,7) = neighbor.mu;
    R(7,7) = -std::sqrt(neighbor.mu / neighbor.rho);

    R(0,8) = neighbor.lambda + 2*neighbor.mu;
    R(1,8) = neighbor.lambda;
    R(2,8) = neighbor.lambda;
    R(6,8) = -std::sqrt((neighbor.lambda + 2 * neighbor.mu) / neighbor.rho);
  }


  if (faceType == FaceType::freeSurface) {
    getTransposedFreeSurfaceGodunovState(local, QgodLocal, QgodNeighbor, R);
  } else {
    Matrix99 chi = Matrix99::Zero();
    chi(0,0) = 1.0;
    chi(1,1) = 1.0;
    chi(2,2) = 1.0;

    const auto godunov = ((R*chi)*R.inverse()).eval();

    // QgodLocal = I - QgodNeighbor
    for (unsigned i = 0; i < godunov.cols(); ++i) {
      for (unsigned j = 0; j < godunov.rows(); ++j) {
        QgodLocal(i,j) = -godunov(j,i);
        QgodNeighbor(i,j) = godunov(j,i);
      }
    }
    for (unsigned idx = 0; idx < 9; ++idx) {
      QgodLocal(idx,idx) += 1.0;
    }
  }
}
#endif // MODEL_COMMON_HPP_
