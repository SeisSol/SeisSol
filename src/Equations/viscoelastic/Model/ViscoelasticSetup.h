// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT in.tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 */

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_VISCOELASTICSETUP_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_VISCOELASTICSETUP_H_

#include "Model/Common.h"
#include "Kernels/Common.h"
#include "Numerical/Transformation.h"

#include <yateto.h>
#include "generated_code/init.h"

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

    template<typename T>
    inline void getTransposedSourceCoefficientTensor( ViscoElasticMaterial const& material, T& sourceMatrix )
    {
      sourceMatrix.setZero();
    
      //       | E_1^T |
      // E^T = |  ...  |
      //       | E_L^T |
      for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
        unsigned offset = 9 + mech * 6;
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
      for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
        for (unsigned i = 0; i < 6; ++i) {
          unsigned idx = 9 + 6*mech + i;
          sourceMatrix(idx, idx) = -material.omega[mech];
        }
      }
    }
    
    template<typename T>
    inline void getTransposedCoefficientMatrix(ViscoElasticMaterial const& i_material, unsigned i_dim, T& AT)
    {
      getTransposedCoefficientMatrix(dynamic_cast<ElasticMaterial const&>(i_material), i_dim, AT);
    
      for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
        getTransposedViscoelasticCoefficientMatrix( i_material.omega[mech],
                                                    i_dim,
                                                    mech,
                                                    AT );
      }
    }
    
    template<>
    inline void getTransposedGodunovState( ViscoElasticMaterial const& local,
                                           ViscoElasticMaterial const& neighbor,
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

    template<>
    inline void initializeSpecificLocalData( ViscoElasticMaterial const& material,
                                             real timeStepWidth,
                                             ViscoElasticLocalData* localData )
    {
      auto sourceMatrix = init::ET::view::create(localData->sourceMatrix);
      getTransposedSourceCoefficientTensor(material, sourceMatrix);
    }
  }
}


#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_VISCOELASTICSETUP_H_

