/******************************************************************************
** Copyright (c) 2014-2015, Intel Corporation                                **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Alexander Heinecke (Intel Corp.)
******************************************************************************/
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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

#ifndef WAVEPROP_KERNEL_TIME_STPP_H_
#define WAVEPROP_KERNEL_TIME_STPP_H_

#include <Equations/poroelastic/Model/datastructures.hpp>
#include <generated_code/kernel.h>
#include <generated_code/init.h>
#include "Common/constants.hpp"
#include "Equations/datastructures.hpp"
#include "Equations/Time.hpp"

#ifdef ACL_DEVICE
#include <device.h>
#endif // ACL_DEVICE

struct GlobalData;
namespace seissol::waveprop::kernel::time {
template<typename Config>
    class Time<Config, std::enable_if_t<Config::MaterialT::Solver == seissol::model::LocalSolver::SpaceTimePredictorPoroelastic>> {
        using RealT = typename Config::RealT;
  protected:
    static void checkGlobalData(GlobalData const* global, size_t alignment) {
        assert( ((uintptr_t)global->stiffnessMatricesTransposed(0)) % alignment == 0 );
        assert( ((uintptr_t)global->stiffnessMatricesTransposed(1)) % alignment == 0 );
        assert( ((uintptr_t)global->stiffnessMatricesTransposed(2)) % alignment == 0 );
    }

    kernel::spaceTimePredictor m_krnlPrototype;
    kernel::projectDerivativeToNodalBoundaryRotated projectDerivativeToNodalBoundaryRotated;


  /*
   *! Offsets of the derivatives.
   *
   * * Offset counting starts at the zeroth derivative with o_derivativesOffset[0]=0; increasing derivatives follow:
   *   1st derivative: o_derivativesOffset[1]
   *   2nd derivative: o_derivativesOffset[2]
   *   ...
   * * Offset are always counted from position zero; for example the sixth derivative will include all jumps over prior derivatives 0 to 5.
   */
  unsigned int m_derivativesOffsets[Config::ConvergenceOrder];

#ifdef ACL_DEVICE
    kernel::gpu_derivative deviceKrnlPrototype;
    kernel::gpu_projectDerivativeToNodalBoundaryRotated deviceDerivativeToNodalBoundaryRotated;
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif

public:
    Time(){
  m_derivativesOffsets[0] = 0;
  for (int order = 0; order < ConvergenceOrder; ++order) {
    if (order > 0) {
      m_derivativesOffsets[order] = tensor::dQ::size(order-1) + m_derivativesOffsets[order-1];
    }
  }
}

void setHostGlobalData(GlobalData const* global) {
  for (int n = 0; n < ConvergenceOrder; ++n) {
    if (n > 0) {
      for (int d = 0; d < 3; ++d) {
        m_krnlPrototype.kDivMTSub(d,n) = init::kDivMTSub::Values[tensor::kDivMTSub::index(d,n)];
      }
    }
    m_krnlPrototype.selectModes(n) = init::selectModes::Values[tensor::selectModes::index(n)];
  }
  for (int k = 0; k < NUMBER_OF_QUANTITIES; k++) {
    m_krnlPrototype.selectQuantity(k) = init::selectQuantity::Values[tensor::selectQuantity::index(k)];
    m_krnlPrototype.selectQuantityG(k) = init::selectQuantityG::Values[tensor::selectQuantityG::index(k)];
  }
  m_krnlPrototype.timeInt = init::timeInt::Values;
  m_krnlPrototype.wHat = init::wHat::Values;
}

void setGlobalData(const CompoundGlobalData& global) {
  setHostGlobalData(global.onHost);

#ifdef ACL_DEVICE
  logError() << "Poroelasticity does not work on GPUs.";
#endif
}

void executeSTP( double                      i_timeStepWidth,
                                         LocalData<Config>&                  data,
                                         RealT                        o_timeIntegrated[tensor::I::size()],
                                         RealT*                       stp )

{
  alignas(PAGESIZE_STACK) RealT stpRhs[tensor::spaceTimePredictorRhs::size()];
  assert( ((uintptr_t)stp) % Alignment  == 0);
  std::fill(std::begin(stpRhs), std::end(stpRhs), 0);
  std::fill(stp, stp + tensor::spaceTimePredictor::size(), 0);
  kernel::spaceTimePredictor krnl = m_krnlPrototype;
 
  //libxsmm can not generate GEMMs with alpha!=1. As a workaround we multiply the 
  //star matrices with dt before we execute the kernel.
  RealT A_values[init::star::size(0)];
  RealT B_values[init::star::size(1)];
  RealT C_values[init::star::size(2)];
  for (size_t i = 0; i < init::star::size(0); i++) {
    A_values[i] = i_timeStepWidth * data.localIntegration.starMatrices[0][i];
    B_values[i] = i_timeStepWidth * data.localIntegration.starMatrices[1][i];
    C_values[i] = i_timeStepWidth * data.localIntegration.starMatrices[2][i];
  }
  krnl.star(0) = A_values;
  krnl.star(1) = B_values;
  krnl.star(2) = C_values;

  //The matrix Zinv depends on the timestep
  //If the timestep is not as expected e.g. when approaching a sync point
  //we have to recalculate it
  if (i_timeStepWidth != data.localIntegration.specific.typicalTimeStepWidth) {
    auto sourceMatrix = init::ET::view::create(data.localIntegration.specific.sourceMatrix);
    RealT ZinvData[NUMBER_OF_QUANTITIES][ConvergenceOrder*ConvergenceOrder];
    model::zInvInitializerForLoop<0, NUMBER_OF_QUANTITIES, decltype(sourceMatrix)>(ZinvData, sourceMatrix, i_timeStepWidth);
    for (size_t i = 0; i < NUMBER_OF_QUANTITIES; i++) {
      krnl.Zinv(i) = ZinvData[i];
    }
  } else {
    for (size_t i = 0; i < NUMBER_OF_QUANTITIES; i++) {
      krnl.Zinv(i) = data.localIntegration.specific.Zinv[i];
    }
  }
  krnl.Gk = data.localIntegration.specific.G[10] * i_timeStepWidth;
  krnl.Gl = data.localIntegration.specific.G[11] * i_timeStepWidth;
  krnl.Gm = data.localIntegration.specific.G[12] * i_timeStepWidth;

  krnl.Q = const_cast<RealT*>(data.dofs);
  krnl.I = o_timeIntegrated;
  krnl.timestep = i_timeStepWidth;
  krnl.spaceTimePredictor = stp;
  krnl.spaceTimePredictorRhs = stpRhs;
  krnl.execute();
}
                                          

void computeAder( double i_timeStepWidth,
                                          LocalData<Config>& data,
                                          LocalTmp<Config>& tmp,
                                          RealT o_timeIntegrated[tensor::I::size()],
                                          RealT* o_timeDerivatives,
                                          bool updateDisplacement)
{
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)data.dofs)              % Alignment == 0 );
  assert( ((uintptr_t)o_timeIntegrated )      % Alignment == 0 );
  assert( ((uintptr_t)o_timeDerivatives)      % Alignment == 0 || o_timeDerivatives == NULL );

  alignas(PAGESIZE_STACK) RealT temporaryBuffer[tensor::spaceTimePredictor::size()];
  RealT* stpBuffer = (o_timeDerivatives != nullptr) ? o_timeDerivatives : temporaryBuffer;
  executeSTP( i_timeStepWidth, data, o_timeIntegrated, stpBuffer );
}

void flopsAder( unsigned int        &o_nonZeroFlops,
                                        unsigned int        &o_hardwareFlops ) {
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops =0;

  o_nonZeroFlops = kernel::spaceTimePredictor::NonZeroFlops;
  o_hardwareFlops = kernel::spaceTimePredictor::HardwareFlops;
  //we multiply the star matrices with dt before we execute the kernel
  o_nonZeroFlops += 3*init::star::size(0);
  o_hardwareFlops += 3*init::star::size(0);
}

unsigned bytesAder()
{
  unsigned reals = 0;
  
  // DOFs load, tDOFs load, tDOFs write
  reals += tensor::Q::size() + 2 * tensor::I::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star>();
  // Zinv
  reals += yateto::computeFamilySize<tensor::Zinv>();
  // G
  reals += 3;
           
  /// \todo incorporate derivatives

  return reals * sizeof(RealT);
}

void computeIntegral( double                            i_expansionPoint,
                                              double                            i_integrationStart,
                                              double                            i_integrationEnd,
                                              const RealT*                       i_timeDerivatives,
                                              RealT                              o_timeIntegrated[tensor::I::size()])
{
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)i_timeDerivatives)  % Alignment == 0 );
  assert( ((uintptr_t)o_timeIntegrated)   % Alignment == 0 );

  // assert that this is a forwared integration in time
  assert( i_integrationStart + (RealT) 1.E-10 > i_expansionPoint   );
  assert( i_integrationEnd                   > i_integrationStart );

  /*
   * compute time integral.
   */
  // compute lengths of integration intervals
  RealT l_deltaTLower = i_integrationStart - i_expansionPoint;
  RealT l_deltaTUpper = i_integrationEnd   - i_expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  RealT l_firstTerm  = (RealT) 1;
  RealT l_secondTerm = (RealT) 1;
  RealT l_factorial  = (RealT) 1;
  
  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = o_timeIntegrated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = i_timeDerivatives + m_derivativesOffsets[i];
  }
 
  // iterate over time derivatives
  for(int der = 0; der < ConvergenceOrder; ++der ) {
    l_firstTerm  *= l_deltaTUpper;
    l_secondTerm *= l_deltaTLower;
    l_factorial  *= (RealT)(der+1);

    intKrnl.power  = l_firstTerm - l_secondTerm;
    intKrnl.power /= l_factorial;

    intKrnl.execute(der);
  }
}

void computeTaylorExpansion( RealT         time,
                                                     RealT         expansionPoint,
                                                     RealT const*  timeDerivatives,
                                                     RealT         timeEvaluated[tensor::Q::size()] ) {
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)timeDerivatives)  % Alignment == 0 );
  assert( ((uintptr_t)timeEvaluated)    % Alignment == 0 );

  // assert that this is a forward evaluation in time
  assert( time >= expansionPoint );

  RealT deltaT = time - expansionPoint;

  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");

  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = timeEvaluated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = timeDerivatives + m_derivativesOffsets[i];
  }
  intKrnl.power = 1.0;
 
  // iterate over time derivatives
  for(int derivative = 0; derivative < ConvergenceOrder; ++derivative) {
    intKrnl.execute(derivative);
    intKrnl.power *= deltaT / RealT(derivative+1);
  }
}

void flopsTaylorExpansion(long long& nonZeroFlops, long long& hardwareFlops) {
  // reset flops
  nonZeroFlops = 0; hardwareFlops = 0;

  // interate over derivatives
  for (unsigned der = 0; der < ConvergenceOrder; ++der) {
    nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(der);
    hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(der);
  }
}

};

}

#endif

