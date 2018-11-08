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
#include <generated_code/sizes.h>

namespace seissol {
  namespace initializers {
    struct DynamicRupture;
  }
}

struct seissol::initializers::DynamicRupture {
  Variable<real*>                                       timeDerivativePlus;
  Variable<real*>                                       timeDerivativeMinus;
  Variable<real[CONVERGENCE_ORDER][seissol::model::godunovState::reals]>   godunov;
  Variable<real[seissol::model::godunovState::reals]>   imposedStatePlus;
  Variable<real[seissol::model::godunovState::reals]>   imposedStateMinus;
  Variable<DRGodunovData>                               godunovData;
  Variable<real[seissol::model::fluxSolver::reals]>     fluxSolverPlus;
  Variable<real[seissol::model::fluxSolver::reals]>     fluxSolverMinus;
  Variable<DRFaceInformation>                           faceInformation;
  Variable<seissol::model::IsotropicWaveSpeeds>         waveSpeedsPlus;
  Variable<seissol::model::IsotropicWaveSpeeds>         waveSpeedsMinus;
  Variable<real[seissol::model::godunovState::rows]>    absoluteSlip;
  
  
  void addTo(LTSTree& tree) {
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(      timeDerivativePlus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(     timeDerivativeMinus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(                 godunov,             mask,     PAGESIZE_HEAP,      seissol::memory::Standard );
    tree.addVar(        imposedStatePlus,             mask,     PAGESIZE_HEAP,      seissol::memory::Standard );
    tree.addVar(       imposedStateMinus,             mask,     PAGESIZE_HEAP,      seissol::memory::Standard );
    tree.addVar(             godunovData,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(          fluxSolverPlus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(         fluxSolverMinus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(         faceInformation,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(          waveSpeedsPlus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(         waveSpeedsMinus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(            absoluteSlip,             mask,         ALIGNMENT,      seissol::memory::Standard );
  }
};
#endif
