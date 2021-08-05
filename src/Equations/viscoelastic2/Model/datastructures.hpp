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

#ifndef MODEL_VISCOELASTIC2_DATASTRUCTURES_H_
#define MODEL_VISCOELASTIC2_DATASTRUCTURES_H_

#include <Model/common_datastructures.hpp>
#include <Equations/elastic/Model/datastructures.hpp>
#include <generated_code/tensor.h>

namespace seissol {
  namespace model {

    struct ViscoElasticMaterial : public ElasticMaterial {
      //! Relaxation frequencies
      double omega[ALLOW_POSSILBE_ZERO_LENGTH_ARRAY(NUMBER_OF_RELAXATION_MECHANISMS)];
      /** Entries of the source matrix (E)
       * theta[0] = -(lambda * Y_lambda + 2.0 * mu * Y_mu)
       * theta[1] = -lambda * Y_lambda
       * theta[2] = -2.0 * mu * Y_mu
       **/
      double theta[ALLOW_POSSILBE_ZERO_LENGTH_ARRAY(NUMBER_OF_RELAXATION_MECHANISMS)][3];
      double Qp;
      double Qs;

      ViscoElasticMaterial() {};
      ViscoElasticMaterial( double* materialValues, int numMaterialValues)
      {
        assert(numMaterialValues == 3 + NUMBER_OF_RELAXATION_MECHANISMS * 4);

        this->rho = materialValues[0];
        this->mu = materialValues[1];
        this->lambda = materialValues[2];

        for (int mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
          this->omega[mech] = materialValues[3 + 4*mech];
          for (unsigned i = 1; i < 4; ++i) {
            this->theta[mech][i-1] = materialValues[3 + 4*mech + i];
          }
        }
        //This constructor is used to initialize a ViscoElasticMaterial
        //from the values in Fortran. Qp and Qs are not part of the 
        //material in Fortran, so we set these to NaN.
        Qp = std::numeric_limits<double>::signaling_NaN();
        Qs = std::numeric_limits<double>::signaling_NaN();
      }

      virtual ~ViscoElasticMaterial() {};

      MaterialType getMaterialType() const override {
        return MaterialType::viscoelastic;
      }
    };
  }
}

#endif
