/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013-2015, SeisSol Group
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
 * Time kernel of SeisSol.
 **/

#include "Kernels/Time.h"

#ifndef NDEBUG
#pragma message "compiling time kernel with assertions"
extern long long libxsmm_num_total_flops;
#endif

#include <cstring>
#include <cassert>
#include <stdint.h>

#include <yateto.h>

#include <Kernels/denseMatrixOps.hpp>
#include <generated_code/init.h>

void seissol::kernels::Time::setGlobalData(GlobalData const* global) {
  assert( ((uintptr_t)global->stiffnessMatricesTransposed(0)) % ALIGNMENT == 0 );
  assert( ((uintptr_t)global->stiffnessMatricesTransposed(1)) % ALIGNMENT == 0 );
  assert( ((uintptr_t)global->stiffnessMatricesTransposed(2)) % ALIGNMENT == 0 );

  m_krnlPrototype.kDivMT = global->stiffnessMatricesTransposed;
  m_krnlPrototype.selectAne = init::selectAne::Values;
  m_krnlPrototype.selectEla = init::selectEla::Values;
}

void seissol::kernels::Time::computeAder( double                      i_timeStepWidth,
                                          LocalData&                  data,
                                          LocalTmp&                   tmp,
                                          real                        o_timeIntegrated[tensor::I::size()],
                                          real*                       o_timeDerivatives )
{
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)data.dofs) % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeIntegrated) % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeDerivatives) % ALIGNMENT == 0 || o_timeDerivatives == NULL );

  /*
   * compute ADER scheme.
   */
  // temporary result
  real temporaryBuffer[2][tensor::dQ::size(0)] __attribute__((aligned(PAGESIZE_STACK)));
  real temporaryBufferExt[2][tensor::dQext::size(1)] __attribute__((aligned(PAGESIZE_STACK)));
  real temporaryBufferAne[2][tensor::dQane::size(0)] __attribute__((aligned(PAGESIZE_STACK)));

  kernel::derivative krnl = m_krnlPrototype;
  kernel::derivativeTaylorExpansion intKrnl;

  krnl.dQ(0) = const_cast<real*>(data.dofs);
  intKrnl.dQ(0) = data.dofs;
  if (o_timeDerivatives != nullptr) {
    streamstore(tensor::dQ::size(0), data.dofs, o_timeDerivatives);
    real* derOut = o_timeDerivatives;
    for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
      derOut += tensor::dQ::size(i-1);
      krnl.dQ(i) = derOut;
      intKrnl.dQ(i) = derOut;
    }
  } else {
    for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
      krnl.dQ(i) = temporaryBuffer[i%2];
      intKrnl.dQ(i) = temporaryBuffer[i%2];
    }
  }

  krnl.dQane(0) = const_cast<real*>(data.dofsAne);
  intKrnl.dQane(0) = const_cast<real*>(data.dofsAne);
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQane(i) = temporaryBufferAne[i%2];
    intKrnl.dQane(i) = temporaryBufferAne[i%2];
    krnl.dQext(i) = temporaryBufferExt[i%2];
  }

  intKrnl.I = o_timeIntegrated;
  intKrnl.Iane = tmp.timeIntegratedAne;

  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    krnl.star(i) = data.localIntegration.starMatrices[i];
  }
  krnl.w = data.localIntegration.specific.w;
  krnl.W = data.localIntegration.specific.W;
  krnl.E = data.localIntegration.specific.E;
  
  // powers in the taylor-series expansion
  intKrnl.power = i_timeStepWidth;

  intKrnl.execute0();
  
  for (unsigned der = 1; der < CONVERGENCE_ORDER; ++der) {
    krnl.execute(der);

    // update scalar for this derivative
    intKrnl.power *= i_timeStepWidth / real(der+1);    
    intKrnl.execute(der);
  }
}

void seissol::kernels::Time::flopsAder( unsigned int        &o_nonZeroFlops,
                                        unsigned int        &o_hardwareFlops ) {  
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops =0;

  // initialization
  o_nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(0);
  o_hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(0);

  // interate over derivatives
  for (unsigned der = 1; der < CONVERGENCE_ORDER; ++der) {
    o_nonZeroFlops  += kernel::derivative::nonZeroFlops(der);
    o_hardwareFlops += kernel::derivative::hardwareFlops(der);

    // update of time integrated DOFs
    o_nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(der);
    o_hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(der);
  }
}

unsigned seissol::kernels::Time::bytesAder()
{
  unsigned reals = 0;
  
  // DOFs load, tDOFs load, tDOFs write
  reals += tensor::Q::size() + tensor::Qane::size() + 2 * tensor::I::size() + 2 * tensor::Iane::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star>()
        + tensor::w::size()
        + tensor::W::size()
        + tensor::E::size();
           
  /// \todo incorporate derivatives

  return reals * sizeof(real);
}

void seissol::kernels::Time::computeIntegral( double                                      i_expansionPoint,
                                              double                                      i_integrationStart,
                                              double                                      i_integrationEnd,
                                              real const*                                 i_timeDerivatives,
                                              real                                        o_timeIntegrated[tensor::I::size()] )
{

  /*
   * assert alignments.
   */
  assert( ((uintptr_t)i_timeDerivatives)  % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeIntegrated)   % ALIGNMENT == 0 );

  // compute lengths of integration intervals
  real l_deltaTLower = i_integrationStart - i_expansionPoint;
  real l_deltaTUpper = i_integrationEnd   - i_expansionPoint;
  
  // initialization of scalars in the taylor series expansion (0th term)
  real l_firstTerm  = static_cast<real>(1.0);
  real l_secondTerm = static_cast<real>(1.0);
  real l_factorial  = static_cast<real>(1.0);
  
  kernel::derivativeTaylorExpansionEla intKrnl;
  intKrnl.I = o_timeIntegrated;
  real const* der = i_timeDerivatives;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = der;
    der += tensor::dQ::size(i);
  }

  // iterate over time derivatives
  for(int der = 0; der < CONVERGENCE_ORDER; ++der ) {
    l_firstTerm  *= l_deltaTUpper;
    l_secondTerm *= l_deltaTLower;
    l_factorial  *= static_cast<real>(der+1);

    intKrnl.power  = l_firstTerm - l_secondTerm;
    intKrnl.power /= l_factorial;

    intKrnl.execute(der);
  }
}

void seissol::kernels::Time::computeTaylorExpansion( real         time,
                                                     real         expansionPoint,
                                                     real const*  timeDerivatives,
                                                     real         timeEvaluated[tensor::Q::size()] ) {
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)timeDerivatives)  % ALIGNMENT == 0 );
  assert( ((uintptr_t)timeEvaluated)    % ALIGNMENT == 0 );

  // assert that this is a forward evaluation in time
  assert( time >= expansionPoint );

  real deltaT = time - expansionPoint;

  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");

  kernel::derivativeTaylorExpansionEla intKrnl;
  intKrnl.I = timeEvaluated;
  real const* der = timeDerivatives;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = der;
    der += tensor::dQ::size(i);
  }
  intKrnl.power = 1.0;
 
  // iterate over time derivatives
  for(int derivative = 0; derivative < CONVERGENCE_ORDER; ++derivative) {
    intKrnl.execute(derivative);
    intKrnl.power *= deltaT / real(derivative+1);
  }
}

void seissol::kernels::Time::flopsTaylorExpansion(long long& nonZeroFlops, long long& hardwareFlops) {
  // reset flops
  nonZeroFlops = 0; hardwareFlops = 0;

  // interate over derivatives
  for (unsigned der = 0; der < CONVERGENCE_ORDER; ++der) {
    nonZeroFlops  += kernel::derivativeTaylorExpansionEla::nonZeroFlops(der);
    hardwareFlops += kernel::derivativeTaylorExpansionEla::hardwareFlops(der);
  }
}
