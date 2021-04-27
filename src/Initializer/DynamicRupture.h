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
 
#ifndef INITIALIZER_DR_H_
#define INITIALIZER_DR_H_

#include <Initializer/typedefs.hpp>
#include <Initializer/tree/LTSTree.hpp>
#include <generated_code/tensor.h>

namespace seissol {
  namespace initializers {
    struct DynamicRupture;
  }
}

#ifndef ACL_DEVICE
#	define MEMKIND_NEIGHBOUR_INTEGRATION seissol::memory::Standard
#	define MEMKIND_Q_INTERPOLATED seissol::memory::Standard
#	define MEMKIND_IMPOSED_STATE seissol::memory::Standard
#else
#	define MEMKIND_NEIGHBOUR_INTEGRATION seissol::memory::DeviceUnifiedMemory
#	define MEMKIND_Q_INTERPOLATED seissol::memory::PinnedMemory
#	define MEMKIND_IMPOSED_STATE seissol::memory::DeviceGlobalMemory
#endif

struct seissol::initializers::DynamicRupture {
  Variable<real*>                                                   timeDerivativePlus;
  Variable<real*>                                                   timeDerivativeMinus;
  Variable<real[tensor::QInterpolated::size()]>                     imposedStatePlus;
  Variable<real[tensor::QInterpolated::size()]>                     imposedStateMinus;
  Variable<DRGodunovData>                                           godunovData;
  Variable<real[tensor::fluxSolver::size()]>                        fluxSolverPlus;
  Variable<real[tensor::fluxSolver::size()]>                        fluxSolverMinus;
  Variable<DRFaceInformation>                                       faceInformation;
  Variable<model::IsotropicWaveSpeeds>                              waveSpeedsPlus;
  Variable<model::IsotropicWaveSpeeds>                              waveSpeedsMinus;

#ifdef ACL_DEVICE
  ScratchpadMemory                        idofsPlusOnDevice;
  ScratchpadMemory                        idofsMinusOnDevice;
  ScratchpadMemory                        QInterpolatedPlusOnDevice;
  ScratchpadMemory                        QInterpolatedMinusOnDevice;

  ScratchpadMemory                        QInterpolatedPlusOnHost;
  ScratchpadMemory                        QInterpolatedMinusOnHost;
  ScratchpadMemory                        imposedStatePlusOnHost;
  ScratchpadMemory                        imposedStateMinusOnHost;
#endif
  
  void addTo(LTSTree& tree) {
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(      timeDerivativePlus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(     timeDerivativeMinus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(        imposedStatePlus,             mask,     PAGESIZE_HEAP,      MEMKIND_IMPOSED_STATE );
    tree.addVar(       imposedStateMinus,             mask,     PAGESIZE_HEAP,      MEMKIND_IMPOSED_STATE );
    tree.addVar(             godunovData,             mask,                 1,      MEMKIND_NEIGHBOUR_INTEGRATION );
    tree.addVar(          fluxSolverPlus,             mask,                 1,      MEMKIND_NEIGHBOUR_INTEGRATION );
    tree.addVar(         fluxSolverMinus,             mask,                 1,      MEMKIND_NEIGHBOUR_INTEGRATION );
    tree.addVar(         faceInformation,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(          waveSpeedsPlus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(         waveSpeedsMinus,             mask,                 1,      seissol::memory::Standard );
#ifdef ACL_DEVICE
    tree.addScratchpadMemory(  idofsPlusOnDevice,              1,      seissol::memory::DeviceGlobalMemory);
    tree.addScratchpadMemory(  idofsMinusOnDevice,             1,      seissol::memory::DeviceGlobalMemory);
    tree.addScratchpadMemory(  QInterpolatedPlusOnDevice,      1,      seissol::memory::DeviceGlobalMemory);
    tree.addScratchpadMemory(  QInterpolatedMinusOnDevice,     1,      seissol::memory::DeviceGlobalMemory);

    tree.addScratchpadMemory(  QInterpolatedPlusOnHost,         1,      seissol::memory::PinnedMemory);
    tree.addScratchpadMemory(  QInterpolatedMinusOnHost,        1,      seissol::memory::PinnedMemory);
    tree.addScratchpadMemory(  imposedStatePlusOnHost,          1,      seissol::memory::PinnedMemory);
    tree.addScratchpadMemory(  imposedStateMinusOnHost,         1,      seissol::memory::PinnedMemory);
#endif
  }
};
#endif
