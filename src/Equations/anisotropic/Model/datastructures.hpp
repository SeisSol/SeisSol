/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 20(0,4), SeisSol Group
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

#ifndef MODEL_ANISOTROPIC_DATASTRUCTURES_H_
#define MODEL_ANISOTROPIC_DATASTRUCTURES_H_

#include <Model/common_datastructures.hpp>
#include <Equations/elastic/Model/datastructures.hpp>
#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <generated_code/init.h>
#include <generated_code/tensor.h>
#include <generated_code/kernel.h>

namespace seissol {
  namespace model {
    struct AnisotropicMaterial : Material {
      double c11;
      double c12;
      double c13;
      double c14;
      double c15;
      double c16;
      double c22;
      double c23;
      double c24;
      double c25;
      double c26;
      double c33;
      double c34;
      double c35;
      double c36;
      double c44;
      double c45;
      double c46;
      double c55;
      double c56;
      double c66;

      AnisotropicMaterial() {};

      AnisotropicMaterial(ElasticMaterial m) {
        rho = m.rho;
        c11 = m.lambda + 2*m.mu;
        c12 = m.lambda;
        c13 = m.lambda;
        c14 = 0;
        c15 = 0;
        c16 = 0; 
        c22 = m.lambda + 2*m.mu;
        c23 = m.lambda;
        c24 = 0;
        c25 = 0;
        c26 = 0;
        c33 = m.lambda + 2*m.mu;
        c34 = 0;
        c35 = 0;
        c36 = 0;
        c44 = m.mu; 
        c45 = 0;
        c46 = 0;
        c55 = m.mu; 
        c56 = 0; 
        c66 = m.mu; 
      }

      virtual ~AnisotropicMaterial() {};

      
      void getFullElasticTensor(real fullTensor[81]) const final {
        auto elasticTensorView = init::mElasticTensor::view::create(fullTensor);
        elasticTensorView.setZero();
        elasticTensorView(0,0,0,0) = c11;
        elasticTensorView(0,0,0,1) = c16;
        elasticTensorView(0,0,0,2) = c15;
        elasticTensorView(0,0,1,0) = c16;
        elasticTensorView(0,0,1,1) = c12;
        elasticTensorView(0,0,1,2) = c14;
        elasticTensorView(0,0,2,0) = c15;
        elasticTensorView(0,0,2,1) = c14;
        elasticTensorView(0,0,2,2) = c13;
        elasticTensorView(0,1,0,0) = c16;
        elasticTensorView(0,1,0,1) = c66;
        elasticTensorView(0,1,0,2) = c56;
        elasticTensorView(0,1,1,0) = c66;
        elasticTensorView(0,1,1,1) = c26;
        elasticTensorView(0,1,1,2) = c46;
        elasticTensorView(0,1,2,0) = c56;
        elasticTensorView(0,1,2,1) = c46;
        elasticTensorView(0,1,2,2) = c36;
        elasticTensorView(0,2,0,0) = c15;
        elasticTensorView(0,2,0,1) = c56;
        elasticTensorView(0,2,0,2) = c55;
        elasticTensorView(0,2,1,0) = c56;
        elasticTensorView(0,2,1,1) = c25;
        elasticTensorView(0,2,1,2) = c45;
        elasticTensorView(0,2,2,0) = c55;
        elasticTensorView(0,2,2,1) = c45;
        elasticTensorView(0,2,2,2) = c35;
        elasticTensorView(1,0,0,0) = c16;
        elasticTensorView(1,0,0,1) = c66;
        elasticTensorView(1,0,0,2) = c56;
        elasticTensorView(1,0,1,0) = c66;
        elasticTensorView(1,0,1,1) = c26;
        elasticTensorView(1,0,1,2) = c46;
        elasticTensorView(1,0,2,0) = c56;
        elasticTensorView(1,0,2,1) = c46;
        elasticTensorView(1,0,2,2) = c36;
        elasticTensorView(1,1,0,0) = c12;
        elasticTensorView(1,1,0,1) = c26;
        elasticTensorView(1,1,0,2) = c25;
        elasticTensorView(1,1,1,0) = c26;
        elasticTensorView(1,1,1,1) = c22;
        elasticTensorView(1,1,1,2) = c24;
        elasticTensorView(1,1,2,0) = c25;
        elasticTensorView(1,1,2,1) = c24;
        elasticTensorView(1,1,2,2) = c23;
        elasticTensorView(1,2,0,0) = c14;
        elasticTensorView(1,2,0,1) = c46;
        elasticTensorView(1,2,0,2) = c45;
        elasticTensorView(1,2,1,0) = c46;
        elasticTensorView(1,2,1,1) = c24;
        elasticTensorView(1,2,1,2) = c44;
        elasticTensorView(1,2,2,0) = c45;
        elasticTensorView(1,2,2,1) = c44;
        elasticTensorView(1,2,2,2) = c34;
        elasticTensorView(2,0,0,0) = c15;
        elasticTensorView(2,0,0,1) = c56;
        elasticTensorView(2,0,0,2) = c55;
        elasticTensorView(2,0,1,0) = c56;
        elasticTensorView(2,0,1,1) = c25;
        elasticTensorView(2,0,1,2) = c45;
        elasticTensorView(2,0,2,0) = c55;
        elasticTensorView(2,0,2,1) = c45;
        elasticTensorView(2,0,2,2) = c35;
        elasticTensorView(2,1,0,0) = c14;
        elasticTensorView(2,1,0,1) = c46;
        elasticTensorView(2,1,0,2) = c45;
        elasticTensorView(2,1,1,0) = c46;
        elasticTensorView(2,1,1,1) = c24;
        elasticTensorView(2,1,1,2) = c44;
        elasticTensorView(2,1,2,0) = c45;
        elasticTensorView(2,1,2,1) = c44;
        elasticTensorView(2,1,2,2) = c34;
        elasticTensorView(2,2,0,0) = c13;
        elasticTensorView(2,2,0,1) = c36;
        elasticTensorView(2,2,0,2) = c35;
        elasticTensorView(2,2,1,0) = c36;
        elasticTensorView(2,2,1,1) = c23;
        elasticTensorView(2,2,1,2) = c34;
        elasticTensorView(2,2,2,0) = c35;
        elasticTensorView(2,2,2,1) = c34;
        elasticTensorView(2,2,2,2) = c33;
      }

      double getMaxWaveSpeed() const final{
#ifdef USE_ANISOTROPIC
        real samplingDirectionsData[seissol::tensor::samplingDirections::Size];
        std::copy_n(init::samplingDirections::Values,
            seissol::tensor::samplingDirections::Size,
            samplingDirectionsData);
        auto samplingDirections = init::samplingDirections::view::create(samplingDirectionsData);

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> saes;

        double maxEv = 0;

        real fullTensor[81];
        getFullElasticTensor(fullTensor);
        kernel::computeChristoffel computeChristoffel;
        computeChristoffel.C = fullTensor;

        for(unsigned j = 0; j < 200; ++j)
        {
          real n[3] = { samplingDirections(j, 0),
                        samplingDirections(j, 1),
                        samplingDirections(j, 2)
          };
          real M[9];
          computeChristoffel.n = n;
          computeChristoffel.christoffel = M;
          computeChristoffel.execute();

          saes.compute(Eigen::Matrix<real, 3, 3>(M).cast<double>());
          auto eigenvalues = saes.eigenvalues();
          for(unsigned i = 0; i < 3; ++i) {
            maxEv = eigenvalues(i) > maxEv ? eigenvalues(i) : maxEv;
          }
        }
        return sqrt(maxEv / rho);
#else
        return getPWaveSpeed();
#endif
      }

      double getPWaveSpeed() const final {
        double muBar = (c44 + c55 + c66) / 3.0;
        double lambdaBar = (c11 + c22 + c33) / 3.0 - 2.0*muBar;
        return std::sqrt((lambdaBar + 2*muBar) / rho);
      }

      double getSWaveSpeed() const final {
        double muBar = (c44 + c55 + c66) / 3.0;
        return std::sqrt(muBar / rho);
      }

      MaterialType getMaterialType() const {
        return MaterialType::anisotropic;
      }
    };


  }
}

#endif
