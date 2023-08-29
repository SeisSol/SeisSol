/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
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
 
#ifndef INITIALIZER_LTS_H_
#define INITIALIZER_LTS_H_

#include <Common/configtensor.hpp>
#include <Initializer/typedefs.hpp>
#include <Initializer/tree/LTSTree.hpp>
#include <generated_code/tensor.h>
#include <Kernels/common.hpp>

#include "Initializer/tree/VariableContainer.hpp"
#include "Model/plasticity.hpp"

#ifndef ACL_DEVICE
#   define MEMKIND_GLOBAL   seissol::memory::HighBandwidth
#   define MEMKIND_TIMEDOFS (Config::ConvergenceOrder <= 7 ? seissol::memory::HighBandwidth : seissol::memory::Standard)
#   define MEMKIND_CONSTANT (Config::ConvergenceOrder <= 4 ? seissol::memory::HighBandwidth : seissol::memory::Standard)
#   define MEMKIND_DOFS     (Config::ConvergenceOrder <= 3 ? seissol::memory::HighBandwidth : seissol::memory::Standard)
# define MEMKIND_UNIFIED  seissol::memory::Standard
#else // ACL_DEVICE
#	define MEMKIND_GLOBAL   seissol::memory::Standard
#	define MEMKIND_CONSTANT seissol::memory::Standard
#	define MEMKIND_DOFS     seissol::memory::DeviceUnifiedMemory
#	define MEMKIND_TIMEDOFS seissol::memory::DeviceUnifiedMemory
# define MEMKIND_UNIFIED  seissol::memory::DeviceUnifiedMemory
#endif // ACL_DEVICE

namespace seissol {
  namespace tensor {
    class Qane;
  }
}

namespace seissol::initializers {
  template<typename Config>
  struct LTS : LTSVariableContainer {
    using ConfigT = Config;
    using RealT = typename Config::RealT;
    using MaterialT = typename Config::MaterialT;

    Variable<RealT[ConfigConstants<Config>::DofsElaSize]>       dofs;
    // size is zero if Qane is not defined
    Variable<RealT[ZeroLengthArrayHandler(ConfigConstants<Config>::DofsAneSize)]> dofsAne;
    Variable<RealT*>                         buffers;
    Variable<RealT*>                         derivatives;
    Variable<CellLocalInformation>          cellInformation;
    Variable<SecondaryCellLocalInformation>          secondaryCellInformation;
    Variable<RealT*[4]>                      faceNeighbors;
    Variable<LocalIntegrationData<Config>>          localIntegration;
    Variable<NeighboringIntegrationData<Config>>    neighboringIntegration;
    Variable<CellMaterialData>              material;
    Variable<MaterialT>                    materialData;
    Variable<seissol::model::PlasticityData<Config>>                plasticity;
    Variable<CellDRMapping[4]>              drMapping;
    Variable<CellBoundaryMapping[4]>        boundaryMapping;
    Variable<RealT[ConfigConstants<Config>::PStrainSize]> pstrain;
    Variable<RealT*[4]>                      faceDisplacements;
    Bucket                                  buffersDerivatives;
    Bucket                                  faceDisplacementsBuffer;

  #ifdef ACL_DEVICE
    Variable<LocalIntegrationData>          localIntegrationOnDevice;
    Variable<NeighboringIntegrationData>    neighIntegrationOnDevice;
    ScratchpadMemory                        integratedDofsScratch;
    ScratchpadMemory                        derivativesScratch;
    ScratchpadMemory                        nodalAvgDisplacements;
  #endif
    
    /// \todo Memkind
    void addTo(LTSTree& tree) override {
      LayerMask plasticityMask;
      if constexpr (Config::Plasticity) {
        plasticityMask = LayerMask(Ghost);
      } else {
        plasticityMask = LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior);
      }

      tree.addVar(                    dofs, LayerMask(Ghost),     PagesizeHeap,      MEMKIND_DOFS );
      if constexpr (ConfigConstants<Config>::DofsAneSize > 0) {
        tree.addVar(                 dofsAne, LayerMask(Ghost),     PagesizeHeap,      MEMKIND_DOFS );
      }
      tree.addVar(                 buffers,      LayerMask(),                 1,      MEMKIND_TIMEDOFS );
      tree.addVar(             derivatives,      LayerMask(),                 1,      MEMKIND_TIMEDOFS );
      tree.addVar(         cellInformation,      LayerMask(),                 1,      MEMKIND_CONSTANT );
      tree.addVar(         secondaryCellInformation,      LayerMask(),                 1,      MEMKIND_CONSTANT );
      tree.addVar(           faceNeighbors, LayerMask(Ghost),                 1,      MEMKIND_TIMEDOFS );
      tree.addVar(        localIntegration, LayerMask(Ghost),                 1,      MEMKIND_CONSTANT );
      tree.addVar(  neighboringIntegration, LayerMask(Ghost),                 1,      MEMKIND_CONSTANT );
      tree.addVar(                material, LayerMask(Ghost),                 1,      seissol::memory::Standard );
      tree.addVar(                materialData,  LayerMask(),                 1,      seissol::memory::Standard );
      tree.addVar(              plasticity,   plasticityMask,                 1,      MEMKIND_UNIFIED );
      tree.addVar(               drMapping, LayerMask(Ghost),                 1,      MEMKIND_CONSTANT );
      tree.addVar(         boundaryMapping, LayerMask(Ghost),                 1,      MEMKIND_CONSTANT );
      tree.addVar(                 pstrain,   plasticityMask,     PagesizeHeap,      MEMKIND_UNIFIED );
      tree.addVar(       faceDisplacements, LayerMask(Ghost),     PagesizeHeap,      seissol::memory::Standard );

      tree.addBucket(buffersDerivatives,                          PagesizeHeap,      MEMKIND_TIMEDOFS );
      tree.addBucket(faceDisplacementsBuffer,                     PagesizeHeap,      MEMKIND_TIMEDOFS );

  #ifdef ACL_DEVICE
      tree.addVar(   localIntegrationOnDevice,   LayerMask(Ghost),  1,      seissol::memory::DeviceGlobalMemory);
      tree.addVar(   neighIntegrationOnDevice,   LayerMask(Ghost),  1,      seissol::memory::DeviceGlobalMemory);
      tree.addScratchpadMemory(  integratedDofsScratch,             1,      seissol::memory::DeviceUnifiedMemory);
      tree.addScratchpadMemory(derivativesScratch,                  1,      seissol::memory::DeviceGlobalMemory);
      tree.addScratchpadMemory(nodalAvgDisplacements,               1,      seissol::memory::DeviceGlobalMemory);
  #endif
    }
  };
}
#endif
