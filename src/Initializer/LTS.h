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

#include "tree/Layer.hpp"
#include "Initializer/typedefs.hpp"
#include "Initializer/tree/LTSTree.hpp"
#include "generated_code/tensor.h"
#include "Kernels/common.hpp"
#include "Parallel/Helper.hpp"

#ifndef ACL_DEVICE
#   define MEMKIND_GLOBAL   AllocationMode::HostOnlyHBM
#if CONVERGENCE_ORDER <= 7
#   define MEMKIND_TIMEDOFS AllocationMode::HostOnlyHBM
#   define MEMKIND_TIMEDOFS_CONSTANT AllocationMode::HostOnlyHBM
#else
#   define MEMKIND_TIMEDOFS AllocationMode::HostOnly
#   define MEMKIND_TIMEDOFS_CONSTANT AllocationMode::HostOnly
#endif
#define MEMKIND_TIMEBUCKET MEMKIND_TIMEDOFS
#if CONVERGENCE_ORDER <= 4
#   define MEMKIND_CONSTANT AllocationMode::HostOnlyHBM
#else
#   define MEMKIND_CONSTANT AllocationMode::HostOnly
#endif
#if CONVERGENCE_DOFS <= 3
#   define MEMKIND_DOFS     AllocationMode::HostOnlyHBM
#else
#   define MEMKIND_DOFS     AllocationMode::HostOnly
#endif
# define MEMKIND_UNIFIED  AllocationMode::HostOnly

#   define MEMKIND_TIMEDOFS_CONSTANT MEMKIND_TIMEDOFS
#   define MEMKIND_CONSTANT_SHARED MEMKIND_CONSTANT
#else // ACL_DEVICE
#	define MEMKIND_GLOBAL   AllocationMode::HostOnly
#	define MEMKIND_CONSTANT AllocationMode::HostOnly
#	define MEMKIND_CONSTANT_SHARED useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit
# define MEMKIND_TIMEDOFS_CONSTANT AllocationMode::HostOnly
#	define MEMKIND_DOFS     useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit
#	define MEMKIND_TIMEDOFS useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit
#	define MEMKIND_TIMEBUCKET useMPIUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplitPinned
# define MEMKIND_UNIFIED  useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit
#endif // ACL_DEVICE

namespace seissol {
  namespace initializer {
    struct LTS;
  }
  namespace tensor {
    class Qane;
  }
}

struct seissol::initializer::LTS {
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
  Variable<real[tensor::QStress::size() + tensor::QEtaModal::size()]> pstrain;
  Variable<real*[4]>                      faceDisplacements;
  Bucket                                  buffersDerivatives;
  Bucket                                  faceDisplacementsBuffer;

  Variable<real*>                         buffersDevice;
  Variable<real*>                         derivativesDevice;
  Variable<real*[4]>                      faceDisplacementsDevice;
  Variable<real*[4]>                      faceNeighborsDevice;
  Variable<CellDRMapping[4]>              drMappingDevice;
  Variable<CellBoundaryMapping[4]>        boundaryMappingDevice;

#ifdef ACL_DEVICE
  ScratchpadMemory                        integratedDofsScratch;
  ScratchpadMemory                        derivativesScratch;
  ScratchpadMemory                        nodalAvgDisplacements;
  ScratchpadMemory                        derivativesExtScratch;
  ScratchpadMemory                        derivativesAneScratch;
  ScratchpadMemory                        idofsAneScratch;
  ScratchpadMemory                        dofsExtScratch;
  ScratchpadMemory                        analyticScratch;
#endif
  
  /// \todo Memkind
  void addTo(LTSTree& tree, bool usePlasticity) {
    LayerMask plasticityMask;
    if (usePlasticity) {
      plasticityMask = LayerMask(Ghost);
    } else {
      plasticityMask = LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior);
    }

    tree.addVar(                    dofs, LayerMask(Ghost),     PAGESIZE_HEAP,      MEMKIND_DOFS );
    if (kernels::size<tensor::Qane>() > 0) {
      tree.addVar(                 dofsAne, LayerMask(Ghost),     PAGESIZE_HEAP,      MEMKIND_DOFS );
    }
    tree.addVar(                 buffers,      LayerMask(),                 1,      MEMKIND_TIMEDOFS_CONSTANT, true );
    tree.addVar(             derivatives,      LayerMask(),                 1,      MEMKIND_TIMEDOFS_CONSTANT, true );
    tree.addVar(         cellInformation,      LayerMask(),                 1,      MEMKIND_CONSTANT, true );
    tree.addVar(           faceNeighbors, LayerMask(Ghost),                 1,      MEMKIND_TIMEDOFS_CONSTANT, true );
    tree.addVar(        localIntegration, LayerMask(Ghost),                 1,      MEMKIND_CONSTANT_SHARED, true );
    tree.addVar(  neighboringIntegration, LayerMask(Ghost),                 1,      MEMKIND_CONSTANT_SHARED, true );
    tree.addVar(                material, LayerMask(Ghost),                 1,      AllocationMode::HostOnly, true );
    tree.addVar(              plasticity,   plasticityMask,                 1,      MEMKIND_UNIFIED, true );
    tree.addVar(               drMapping, LayerMask(Ghost),                 1,      MEMKIND_CONSTANT, true );
    tree.addVar(         boundaryMapping, LayerMask(Ghost),                 1,      MEMKIND_CONSTANT, true );
    tree.addVar(                 pstrain,   plasticityMask,     PAGESIZE_HEAP,      MEMKIND_UNIFIED );
    tree.addVar(       faceDisplacements, LayerMask(Ghost),     PAGESIZE_HEAP,      AllocationMode::HostOnly, true );

    tree.addBucket(buffersDerivatives,                          PAGESIZE_HEAP,      MEMKIND_TIMEBUCKET );
    tree.addBucket(faceDisplacementsBuffer,                     PAGESIZE_HEAP,      MEMKIND_TIMEDOFS );

    tree.addVar(   buffersDevice, LayerMask(),     1,      AllocationMode::HostOnly, true );
    tree.addVar(   derivativesDevice, LayerMask(),     1,      AllocationMode::HostOnly, true );
    tree.addVar(   faceDisplacementsDevice, LayerMask(Ghost),     1,      AllocationMode::HostOnly, true );
    tree.addVar(   faceNeighborsDevice, LayerMask(Ghost),     1,      AllocationMode::HostOnly, true );
    tree.addVar(   drMappingDevice, LayerMask(Ghost),     1,      AllocationMode::HostOnly, true );
    tree.addVar(   boundaryMappingDevice, LayerMask(Ghost),     1,      AllocationMode::HostOnly, true );

#ifdef ACL_DEVICE
    tree.addScratchpadMemory(derivativesExtScratch,               1,      AllocationMode::DeviceOnly);
    tree.addScratchpadMemory(derivativesAneScratch,               1,      AllocationMode::DeviceOnly);
    tree.addScratchpadMemory(idofsAneScratch,               1,      AllocationMode::DeviceOnly);
    tree.addScratchpadMemory(dofsExtScratch,               1,      AllocationMode::DeviceOnly);
    tree.addScratchpadMemory(  integratedDofsScratch,             1,      AllocationMode::HostDeviceSplit);
    tree.addScratchpadMemory(derivativesScratch,                  1,      AllocationMode::DeviceOnly);
    tree.addScratchpadMemory(nodalAvgDisplacements,               1,      AllocationMode::DeviceOnly);
    tree.addScratchpadMemory(analyticScratch,               1,      AllocationMode::HostDevicePinned);
#endif
  }
};
#endif
