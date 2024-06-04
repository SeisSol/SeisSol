/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2013-2014, SeisSol Group
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
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Volume kernel of SeisSol.
 **/

#ifndef VOLUME_H_
#define VOLUME_H_

#include <Initializer/LTS.h>
#include <Initializer/Parameters/ModelParameters.h>
#include <Initializer/typedefs.hpp>
#include <Kernels/Interface.hpp>
#include <Kernels/LocalBase.h>
#include <Kernels/common.hpp>
#include <cassert>
#include <generated_code/tensor.h>

namespace seissol {
  class SeisSol;
} // namespace seissol
namespace seissol::kernels {
  class Local;
} // namespace seissol::kernels

class seissol::kernels::Local : public LocalBase {
  protected:
  seissol::initializer::parameters::DamagedElasticParameters* m_damagedElasticParameters;

  public:
    void setHostGlobalData(GlobalData const* global);
    void setGlobalData(const CompoundGlobalData& global);

    void computeIntegral(real i_timeIntegratedDegreesOfFreedom[tensor::I::size()],
                         LocalData& data,
                         LocalTmp& tmp,
                         const CellMaterialData* materialData,
                         CellBoundaryMapping const (*cellBoundaryMapping)[4],
                         double time,
                         double timeStepWidth);

    void computeBatchedIntegral(ConditionalPointersToRealsTable& dataTable,
                                ConditionalMaterialTable& materialTable,
                                ConditionalIndicesTable& indicesTable,
                                kernels::LocalData::Loader& loader,
                                LocalTmp& tmp,
                                double timeStepWidth);

    void evaluateBatchedTimeDependentBc(ConditionalPointersToRealsTable& dataTable,
                                        ConditionalIndicesTable& indicesTable,
                                        kernels::LocalData::Loader& loader,
                                        double time,
                                        double timeStepWidth);

    void flopsIntegral(FaceType const i_faceTypes[4],
                       unsigned int &o_nonZeroFlops,
                       unsigned int &o_hardwareFlops );
                        
    unsigned bytesIntegral();

    void setDamagedElasticParameters(seissol::initializer::parameters::DamagedElasticParameters* damagedElasticParameters) {
      m_damagedElasticParameters = damagedElasticParameters;
    };
    #ifdef USE_DAMAGEDELASTIC
      void computeNonLinearRusanovFlux(const CellMaterialData* materialData,
                                   unsigned int l_cell,
                                   unsigned int side,
                                   const double* timeWeights,
                                   const real* qIPlus,
                                   const real* qIMinus,
                                   real* rusanovFluxP,
                                   const LocalIntegrationData* localIntegration);
      void computeNonLinearIntegralCorrection(
      const CellLocalInformation* cellInformation,
      unsigned int l_cell,
      real** derivatives,
      real* (*faceNeighbors)[4],
      const CellMaterialData* materialData,
      const LocalIntegrationData* localIntegration,
      const NeighborData& data,
      const CellDRMapping (*drMapping)[4],
      kernel::nonlinearSurfaceIntegral& m_nonlSurfIntPrototype,
      double timeStepSize,
      const kernel::nonlEvaluateAndRotateQAtInterpolationPoints& m_nonlinearInterpolation,
      seissol::SeisSol& seissolInstance,
      initializer::LTS* m_lts,
      double subTimeStart=0.0);                             
    #endif
};

#endif

