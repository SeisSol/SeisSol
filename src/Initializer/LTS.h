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

#include <Initializer/typedefs.hpp>
#include <Initializer/tree/LTSTree.hpp>
#include <generated_code/tensor.h>
#include <Kernels/common.hpp>

#if CONVERGENCE_ORDER < 2 || CONVERGENCE_ORDER > 8
#error Preprocessor flag CONVERGENCE_ORDER is not in {2, 3, 4, 5, 6, 7, 8}.
#endif

#ifndef ACL_DEVICE
#   define MEMKIND_GLOBAL   seissol::memory::HighBandwidth
#if CONVERGENCE_ORDER <= 7
#   define MEMKIND_TIMEDOFS seissol::memory::HighBandwidth
#else
#   define MEMKIND_TIMEDOFS seissol::memory::Standard
#endif
#if CONVERGENCE_ORDER <= 4
#   define MEMKIND_CONSTANT seissol::memory::HighBandwidth
#else
#   define MEMKIND_CONSTANT seissol::memory::Standard
#endif
#if CONVERGENCE_DOFS <= 3
#   define MEMKIND_DOFS     seissol::memory::HighBandwidth
#else
#   define MEMKIND_DOFS     seissol::memory::Standard
#endif
# define MEMKIND_UNIFIED  seissol::memory::Standard
#else // ACL_DEVICE
#	define MEMKIND_GLOBAL   seissol::memory::Standard
#	define MEMKIND_CONSTANT seissol::memory::Standard
#	define MEMKIND_DOFS     seissol::memory::DeviceUnifiedMemory
#	define MEMKIND_TIMEDOFS seissol::memory::DeviceUnifiedMemory
# define MEMKIND_UNIFIED  seissol::memory::DeviceUnifiedMemory
#endif // ACL_DEVICE

namespace seissol {
  namespace initializers {
    struct LTS;
  }
  namespace tensor {
    class Qane;
  }
}

struct seissol::initializers::LTS {
  Variable<real[tensor::Q::size()]>       dofs;
  // size is zero if Qane is not defined
  Variable<real[ALLOW_POSSILBE_ZERO_LENGTH_ARRAY(kernels::size<tensor::Qane>())]> dofsAne;
  Variable<real*>                         buffers;
  Variable<real*>                         derivatives;
  Variable<CellLocalInformation>          cellInformation;
  Variable<real*[4]>                      faceNeighbors;
  Variable<LocalIntegrationData>          localIntegration;
  Variable<NeighboringIntegrationData>    neighboringIntegration;
  Variable<CellMaterialData>              material;
  Variable<PlasticityData>                plasticity;
  Variable<CellDRMapping[4]>              drMapping;
  Variable<CellBoundaryMapping[4]>        boundaryMapping;
  Variable<real[7]>                       pstrain;
  Variable<real*>                         displacements;
  Bucket                                  buffersDerivatives;
  Bucket                                  displacementsBuffer;

#ifdef ACL_DEVICE
  Variable<LocalIntegrationData>          localIntegrationOnDevice;
  Variable<NeighboringIntegrationData>    neighIntegrationOnDevice;
  ScratchpadMemory                        idofsScratch;
  ScratchpadMemory                        derivativesScratch;
#endif
  
  /// \todo Memkind
  void addTo(LTSTree& tree) {
#ifdef USE_PLASTICITY
    LayerMask plasticityMask = LayerMask(Ghost);
#else
    LayerMask plasticityMask = LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior);
#endif
    tree.addVar(                    dofs, LayerMask(Ghost),     PAGESIZE_HEAP,      MEMKIND_DOFS );
    if (kernels::size<tensor::Qane>() > 0) {
      tree.addVar(                 dofsAne, LayerMask(Ghost),     PAGESIZE_HEAP,      MEMKIND_DOFS );
    }
    tree.addVar(                 buffers,      LayerMask(),                 1,      MEMKIND_TIMEDOFS );
    tree.addVar(             derivatives,      LayerMask(),                 1,      MEMKIND_TIMEDOFS );
    tree.addVar(         cellInformation,      LayerMask(),                 1,      MEMKIND_CONSTANT );
    tree.addVar(           faceNeighbors, LayerMask(Ghost),                 1,      MEMKIND_TIMEDOFS );
    tree.addVar(        localIntegration, LayerMask(Ghost),                 1,      MEMKIND_CONSTANT );
    tree.addVar(  neighboringIntegration, LayerMask(Ghost),                 1,      MEMKIND_CONSTANT );
    tree.addVar(                material, LayerMask(Ghost),                 1,      seissol::memory::Standard );
    tree.addVar(              plasticity,   plasticityMask,                 1,      MEMKIND_UNIFIED );
    tree.addVar(               drMapping, LayerMask(Ghost),                 1,      MEMKIND_CONSTANT );
    tree.addVar(         boundaryMapping, LayerMask(Ghost),                 1,      MEMKIND_CONSTANT );
    tree.addVar(                 pstrain,   plasticityMask,     PAGESIZE_HEAP,      MEMKIND_UNIFIED );
    tree.addVar(           displacements, LayerMask(Ghost),     PAGESIZE_HEAP,      seissol::memory::Standard );
    
    tree.addBucket(buffersDerivatives,                          PAGESIZE_HEAP,      MEMKIND_TIMEDOFS );
    tree.addBucket(displacementsBuffer,                         PAGESIZE_HEAP,      MEMKIND_TIMEDOFS );

#ifdef ACL_DEVICE
    tree.addVar(   localIntegrationOnDevice,   LayerMask(Ghost),  1,      seissol::memory::DeviceGlobalMemory );
    tree.addVar(   neighIntegrationOnDevice,   LayerMask(Ghost),  1,      seissol::memory::DeviceGlobalMemory );
    tree.addScratchpadMemory(  idofsScratch,                      1,      seissol::memory::DeviceGlobalMemory);
    tree.addScratchpadMemory(derivativesScratch,                  1,      seissol::memory::DeviceGlobalMemory);
#endif
  }
};
#endif
