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

      
      void getFullElasticTensor(real fullTensor[81]) {
        fullTensor[0]  = static_cast<real>(c11);
        fullTensor[1]  = static_cast<real>(c16);
        fullTensor[2]  = static_cast<real>(c15);
        fullTensor[3]  = static_cast<real>(c16);
        fullTensor[4]  = static_cast<real>(c12);
        fullTensor[5]  = static_cast<real>(c14);
        fullTensor[6]  = static_cast<real>(c15);
        fullTensor[7]  = static_cast<real>(c14);
        fullTensor[8]  = static_cast<real>(c13);
        fullTensor[9]  = static_cast<real>(c16);
        fullTensor[10] = static_cast<real>(c66);
        fullTensor[11] = static_cast<real>(c56);
        fullTensor[12] = static_cast<real>(c66);
        fullTensor[13] = static_cast<real>(c26);
        fullTensor[14] = static_cast<real>(c46);
        fullTensor[15] = static_cast<real>(c56);
        fullTensor[16] = static_cast<real>(c46);
        fullTensor[17] = static_cast<real>(c36);
        fullTensor[18] = static_cast<real>(c15);
        fullTensor[19] = static_cast<real>(c56);
        fullTensor[20] = static_cast<real>(c55);
        fullTensor[21] = static_cast<real>(c56);
        fullTensor[22] = static_cast<real>(c25);
        fullTensor[23] = static_cast<real>(c35);
        fullTensor[24] = static_cast<real>(c55);
        fullTensor[25] = static_cast<real>(c45);
        fullTensor[26] = static_cast<real>(c35);
        fullTensor[27] = static_cast<real>(c16);
        fullTensor[28] = static_cast<real>(c66);
        fullTensor[29] = static_cast<real>(c56);
        fullTensor[30] = static_cast<real>(c66);
        fullTensor[31] = static_cast<real>(c26);
        fullTensor[32] = static_cast<real>(c46);
        fullTensor[33] = static_cast<real>(c56);
        fullTensor[34] = static_cast<real>(c46);
        fullTensor[35] = static_cast<real>(c36);
        fullTensor[36] = static_cast<real>(c12);
        fullTensor[37] = static_cast<real>(c26);
        fullTensor[38] = static_cast<real>(c25);
        fullTensor[39] = static_cast<real>(c26);
        fullTensor[40] = static_cast<real>(c22);
        fullTensor[41] = static_cast<real>(c24);
        fullTensor[42] = static_cast<real>(c25);
        fullTensor[43] = static_cast<real>(c24);
        fullTensor[44] = static_cast<real>(c23);
        fullTensor[45] = static_cast<real>(c14);
        fullTensor[46] = static_cast<real>(c46);
        fullTensor[47] = static_cast<real>(c45);
        fullTensor[48] = static_cast<real>(c46);
        fullTensor[49] = static_cast<real>(c24);
        fullTensor[50] = static_cast<real>(c44);
        fullTensor[51] = static_cast<real>(c45);
        fullTensor[52] = static_cast<real>(c44);
        fullTensor[53] = static_cast<real>(c34);
        fullTensor[54] = static_cast<real>(c15);
        fullTensor[55] = static_cast<real>(c56);
        fullTensor[56] = static_cast<real>(c55);
        fullTensor[57] = static_cast<real>(c56);
        fullTensor[58] = static_cast<real>(c25);
        fullTensor[59] = static_cast<real>(c45);
        fullTensor[60] = static_cast<real>(c55);
        fullTensor[61] = static_cast<real>(c45);
        fullTensor[62] = static_cast<real>(c35);
        fullTensor[63] = static_cast<real>(c14);
        fullTensor[64] = static_cast<real>(c46);
        fullTensor[65] = static_cast<real>(c45);
        fullTensor[66] = static_cast<real>(c46);
        fullTensor[67] = static_cast<real>(c24);
        fullTensor[68] = static_cast<real>(c44);
        fullTensor[69] = static_cast<real>(c45);
        fullTensor[70] = static_cast<real>(c44);
        fullTensor[71] = static_cast<real>(c34);
        fullTensor[72] = static_cast<real>(c14);
        fullTensor[73] = static_cast<real>(c36);
        fullTensor[74] = static_cast<real>(c35);
        fullTensor[75] = static_cast<real>(c36);
        fullTensor[76] = static_cast<real>(c23);
        fullTensor[77] = static_cast<real>(c34);
        fullTensor[78] = static_cast<real>(c35);
        fullTensor[79] = static_cast<real>(c34);
        fullTensor[80] = static_cast<real>(c33);
      }

      double getMaxWaveSpeed() {
#ifdef USE_ANISOTROPIC
        double samplingDirectionsData[seissol::tensor::samplingDirections::Size];
        std::copy_n(init::samplingDirections::Values,
            seissol::tensor::samplingDirections::Size,
            samplingDirectionsData);
        auto samplingDirections = init::samplingDirections::view::create(samplingDirectionsData);

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> saes;

        double maxEv = 0;

        double fullTensor[81];
        getFullElasticTensor(fullTensor);
        kernel::computeChristoffel computeChristoffel;
        computeChristoffel.C = fullTensor;

        for(unsigned j = 0; j < 200; ++j)
        {
          double n[3] = { samplingDirections(j, 0),
                        samplingDirections(j, 1),
                        samplingDirections(j, 2)
          };
          double M[9];
          computeChristoffel.n = n;
          computeChristoffel.christoffel = M;
          computeChristoffel.execute();

          saes.compute(Eigen::Matrix<double, 3, 3>(M));
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

      double getPWaveSpeed() {
        double muBar = (c44 + c55 + c66) / 3.0;
        double lambdaBar = (c11 + c22 + c33) / 3.0 - 2.0*muBar;
        return std::sqrt((lambdaBar + 2*muBar) / rho);
      }

      double getSWaveSpeed() {
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
