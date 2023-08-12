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

#ifndef WAVEPROP_KERNEL_TIME_CK_H_
#define WAVEPROP_KERNEL_TIME_CK_H_

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
    class Time<Config, std::enable_if_t<Config::MaterialT::Solver == seissol::model::LocalSolver::CauchyKovalevski>> {
        using RealT = typename Config::RealT;
  protected:
    static void checkGlobalData(GlobalData const* global, size_t alignment) {
        assert( ((uintptr_t)global->stiffnessMatricesTransposed(0)) % alignment == 0 );
        assert( ((uintptr_t)global->stiffnessMatricesTransposed(1)) % alignment == 0 );
        assert( ((uintptr_t)global->stiffnessMatricesTransposed(2)) % alignment == 0 );
    }

    // TODO(someone): remove or split into STP and CK proper
#ifdef USE_STP
    kernel::spaceTimePredictor m_krnlPrototype;
#else    
    kernel::derivative m_krnlPrototype;
#endif
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
    /**
     * Constructor, which initializes the time kernel.
     **/
    Time() {
        m_derivativesOffsets[0] = 0;
        for (int order = 0; order < ConvergenceOrder; ++order) {
            if (order > 0) {
            m_derivativesOffsets[order] = tensor::dQ::size(order-1) + m_derivativesOffsets[order-1];
            }
        }
    }

    void setHostGlobalData(GlobalData const* global) {
#ifdef USE_STP
  //Note: We could use the space time predictor for elasticity.
  //This is not tested and experimental
  for (int n = 0; n < ConvergenceOrder; ++n) {
    if (n > 0) {
      for (int d = 0; d < 3; ++d) {
        m_krnlPrototype.kDivMTSub(d,n) = init::kDivMTSub::Values[tensor::kDivMTSub::index(d,n)];
      }
    }
    m_krnlPrototype.selectModes(n) = init::selectModes::Values[tensor::selectModes::index(n)];
  }
  m_krnlPrototype.Zinv = init::Zinv::Values;
  m_krnlPrototype.timeInt = init::timeInt::Values;
  m_krnlPrototype.wHat = init::wHat::Values;
#else //USE_STP
  checkGlobalData(global, Alignment);

  m_krnlPrototype.kDivMT = global->stiffnessMatricesTransposed;

  projectDerivativeToNodalBoundaryRotated.V3mTo2nFace = global->V3mTo2nFace;

#endif //USE_STP
}

void setGlobalData(const CompoundGlobalData& global) {
  setHostGlobalData(global.onHost);

#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);
  const auto deviceAlignment = device.api->getGlobMemAlignment();
  checkGlobalData(global.onDevice, deviceAlignment);
  deviceKrnlPrototype.kDivMT = global.onDevice->stiffnessMatricesTransposed;
  deviceDerivativeToNodalBoundaryRotated.V3mTo2nFace = global.onDevice->V3mTo2nFace;

#endif
}

void computeAder(double i_timeStepWidth,
                                         LocalData<Config>& data,
                                         LocalTmp<Config>& tmp,
                                         RealT o_timeIntegrated[tensor::I::size()],
                                         RealT* o_timeDerivatives,
                                         bool updateDisplacement) {

  assert(reinterpret_cast<uintptr_t>(data.dofs) % Alignment == 0 );
  assert(reinterpret_cast<uintptr_t>(o_timeIntegrated) % Alignment == 0 );
  assert(o_timeDerivatives == nullptr || reinterpret_cast<uintptr_t>(o_timeDerivatives) % Alignment == 0);

  // Only a small fraction of cells has the gravitational free surface boundary condition
  updateDisplacement &= std::any_of(std::begin(data.cellInformation.faceTypes),
                                    std::end(data.cellInformation.faceTypes),
                                    [](const FaceType f) {
                                      return f == FaceType::freeSurfaceGravity;
                                    });

#ifdef USE_STP
  //Note: We could use the space time predictor for elasticity.
  //This is not tested and experimental
  alignas(PAGESIZE_STACK) RealT stpRhs[tensor::spaceTimePredictor::size()];
  alignas(PAGESIZE_STACK) RealT stp[tensor::spaceTimePredictor::size()]{};
  kernel::spaceTimePredictor krnl = m_krnlPrototype;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    krnl.star(i) = data.localIntegration.starMatrices[i];
  }
  krnl.Q = const_cast<RealT*>(data.dofs);
  krnl.I = o_timeIntegrated;
  krnl.timestep = i_timeStepWidth;
  krnl.spaceTimePredictor = stp;
  krnl.spaceTimePredictorRhs = stpRhs;
  krnl.execute();
#else //USE_STP
  alignas(PAGESIZE_STACK) RealT temporaryBuffer[yateto::computeFamilySize<tensor::dQ>()];
  auto* derivativesBuffer = (o_timeDerivatives != nullptr) ? o_timeDerivatives : temporaryBuffer;

  kernel::derivative krnl = m_krnlPrototype;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    krnl.star(i) = data.localIntegration.starMatrices[i];
  }

  // Optional source term
  set_ET(krnl, get_ptr_sourceMatrix(data.localIntegration.specific));

  krnl.dQ(0) = const_cast<RealT*>(data.dofs);
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQ(i) = derivativesBuffer + m_derivativesOffsets[i];
  }

  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = o_timeIntegrated;
  intKrnl.dQ(0) = data.dofs;
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = derivativesBuffer + m_derivativesOffsets[i];
  }
  // powers in the taylor-series expansion
  intKrnl.power = i_timeStepWidth;
  intKrnl.execute0();

  if (updateDisplacement) {
    // First derivative if needed later in kernel
    std::copy_n(data.dofs, tensor::dQ::size(0), derivativesBuffer);
  } else if (o_timeDerivatives != nullptr) {
    // First derivative is not needed here but later
    // Hence stream it out
    streamstore(tensor::dQ::size(0), data.dofs, derivativesBuffer);
  }

  for (unsigned der = 1; der < ConvergenceOrder; ++der) {
    krnl.execute(der);

    // update scalar for this derivative
    intKrnl.power *= i_timeStepWidth / RealT(der+1);    
    intKrnl.execute(der);
  }

  // Do not compute it like this if at interface
  // Compute integrated displacement over time step if needed.
  if (updateDisplacement) {
    auto& bc = tmp.gravitationalFreeSurfaceBc;
    for (unsigned face = 0; face < 4; ++face) {
      if (data.faceDisplacements[face] != nullptr
          && data.cellInformation.faceTypes[face] == FaceType::freeSurfaceGravity) {
        bc.evaluate(
            face,
            projectDerivativeToNodalBoundaryRotated,
            data.boundaryMapping[face],
            data.faceDisplacements[face],
            tmp.nodalAvgDisplacements[face].data(),
            *this,
            derivativesBuffer,
            i_timeStepWidth,
            data.materialData,
            data.cellInformation.faceTypes[face]
        );
      }
    }
  }
#endif //USE_STP
}

void computeBatchedAder(double i_timeStepWidth,
                                                LocalTmp<Config>& tmp,
                                                ConditionalPointersToRealsTable &dataTable,
                                                ConditionalMaterialTable &materialTable,
                                                bool updateDisplacement) {
#ifdef ACL_DEVICE
  kernel::gpu_derivative derivativesKrnl = deviceKrnlPrototype;
  kernel::gpu_derivativeTaylorExpansion intKrnl;

  ConditionalKey timeVolumeKernelKey(KernelNames::Time || KernelNames::Volume);
  if(dataTable.find(timeVolumeKernelKey) != dataTable.end()) {
    auto &entry = dataTable[timeVolumeKernelKey];

    const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
    derivativesKrnl.numElements = numElements;
    intKrnl.numElements = numElements;

    intKrnl.I = (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr();

    unsigned starOffset = 0;
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      derivativesKrnl.star(i) = const_cast<const RealT **>((entry.get(inner_keys::Wp::Id::Star))->getDeviceDataPtr());
      derivativesKrnl.extraOffset_star(i) = starOffset;
      starOffset += tensor::star::size(i);
    }

    unsigned derivativesOffset = 0;
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
      derivativesKrnl.dQ(i) = (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr();
      intKrnl.dQ(i) = const_cast<const RealT **>((entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr());

      derivativesKrnl.extraOffset_dQ(i) = derivativesOffset;
      intKrnl.extraOffset_dQ(i) = derivativesOffset;

      derivativesOffset += tensor::dQ::size(i);
    }

    // stream dofs to the zero derivative
    device.algorithms.streamBatchedData((entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr(),
                                        (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr(),
                                        tensor::Q::Size,
                                        derivativesKrnl.numElements,
                                        device.api->getDefaultStream());

    const auto maxTmpMem = yateto::getMaxTmpMemRequired(intKrnl, derivativesKrnl);
    RealT* tmpMem = reinterpret_cast<RealT*>(device.api->getStackMemory(maxTmpMem * numElements));

    intKrnl.power = i_timeStepWidth;
    intKrnl.linearAllocator.initialize(tmpMem);
    intKrnl.streamPtr = device.api->getDefaultStream();
    intKrnl.execute0();

    for (unsigned Der = 1; Der < ConvergenceOrder; ++Der) {
      derivativesKrnl.linearAllocator.initialize(tmpMem);
      derivativesKrnl.streamPtr = device.api->getDefaultStream();
      derivativesKrnl.execute(Der);

      // update scalar for this derivative
      intKrnl.power *= i_timeStepWidth / RealT(Der + 1);
      intKrnl.linearAllocator.initialize(tmpMem);
      intKrnl.streamPtr = device.api->getDefaultStream();
      intKrnl.execute(Der);
    }
    device.api->popStackMemory();
  }

  if (updateDisplacement) {
    auto& bc = tmp.gravitationalFreeSurfaceBc;
    for (unsigned face = 0; face < 4; ++face) {
      bc.evaluateOnDevice(face,
                          deviceDerivativeToNodalBoundaryRotated,
                          *this,
                          dataTable,
                          materialTable,
                          i_timeStepWidth,
                          device);
    }
  }
#else
  assert(false && "no implementation provided");
#endif
}

void flopsAder( unsigned int        &o_nonZeroFlops,
                                        unsigned int        &o_hardwareFlops ) {
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops =0;

  // initialization
  o_nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(0);
  o_hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(0);

  // interate over derivatives
  for( unsigned l_derivative = 1; l_derivative < ConvergenceOrder; l_derivative++ ) {
    o_nonZeroFlops  += kernel::derivative::nonZeroFlops(l_derivative);
    o_hardwareFlops += kernel::derivative::hardwareFlops(l_derivative);

    // update of time integrated DOFs
    o_nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(l_derivative);
    o_hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(l_derivative);
  }

}

unsigned bytesAder()
{
  unsigned reals = 0;
  
  // DOFs load, tDOFs load, tDOFs write
  reals += tensor::Q::size() + 2 * tensor::I::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star>();
           
  /// \todo incorporate derivatives

  return reals * sizeof(RealT);
}

void computeIntegral( double                            i_expansionPoint,
                                              double                            i_integrationStart,
                                              double                            i_integrationEnd,
                                              const RealT*                       i_timeDerivatives,
                                              RealT                              o_timeIntegrated[tensor::I::size()] )
{
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)i_timeDerivatives)  % Alignment == 0 );
  assert( ((uintptr_t)o_timeIntegrated)   % Alignment == 0 );

  // assert that this is a forwared integration in time
  assert( i_integrationStart + static_cast<RealT>(1.E-10) > i_expansionPoint   );
  assert( i_integrationEnd                   > i_integrationStart );

  /*
   * compute time integral.
   */
  // compute lengths of integration intervals
  RealT l_deltaTLower = i_integrationStart - i_expansionPoint;
  RealT l_deltaTUpper = i_integrationEnd   - i_expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  RealT l_firstTerm  = static_cast<RealT>(1);
  RealT l_secondTerm = static_cast<RealT>(1);
  RealT l_factorial  = static_cast<RealT>(1);
  
  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = o_timeIntegrated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = i_timeDerivatives + m_derivativesOffsets[i];
  }
 
  // iterate over time derivatives
  for(int der = 0; der < ConvergenceOrder; ++der ) {
    l_firstTerm  *= l_deltaTUpper;
    l_secondTerm *= l_deltaTLower;
    l_factorial  *= static_cast<RealT>(der+1);

    intKrnl.power  = l_firstTerm - l_secondTerm;
    intKrnl.power /= l_factorial;

    intKrnl.execute(der);
  }
}

void computeBatchedIntegral(double i_expansionPoint,
                                                    double i_integrationStart,
                                                    double i_integrationEnd,
                                                    const RealT** i_timeDerivatives,
                                                    RealT ** o_timeIntegratedDofs,
                                                    unsigned numElements) {
#ifdef ACL_DEVICE
  // assert that this is a forwared integration in time
  assert( i_integrationStart + static_cast<RealT>(1.E-10) > i_expansionPoint   );
  assert( i_integrationEnd                   > i_integrationStart );

  /*
   * compute time integral.
   */
  // compute lengths of integration intervals
  RealT deltaTLower = i_integrationStart - i_expansionPoint;
  RealT deltaTUpper = i_integrationEnd - i_expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  RealT firstTerm  = static_cast<RealT>(1);
  RealT secondTerm = static_cast<RealT>(1);
  RealT factorial  = static_cast<RealT>(1);

  kernel::gpu_derivativeTaylorExpansion intKrnl;
  intKrnl.numElements = numElements;
  RealT* tmpMem = reinterpret_cast<RealT*>(device.api->getStackMemory(intKrnl.TmpMaxMemRequiredInBytes * numElements));

  intKrnl.I = o_timeIntegratedDofs;

  unsigned derivativesOffset = 0;
  for (size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = i_timeDerivatives;
    intKrnl.extraOffset_dQ(i) = derivativesOffset;
    derivativesOffset += tensor::dQ::size(i);
  }

  // iterate over time derivatives
  for(int der = 0; der < ConvergenceOrder; ++der) {
    firstTerm *= deltaTUpper;
    secondTerm *= deltaTLower;
    factorial *= static_cast<RealT>(der + 1);

    intKrnl.power = firstTerm - secondTerm;
    intKrnl.power /= factorial;
    intKrnl.linearAllocator.initialize(tmpMem);
    intKrnl.streamPtr = device.api->getDefaultStream();
    intKrnl.execute(der);
  }
  device.api->popStackMemory();
#else
  assert(false && "no implementation provided");
#endif
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

void computeBatchedTaylorExpansion(RealT time,
                                                           RealT expansionPoint,
                                                           RealT** timeDerivatives,
                                                           RealT** timeEvaluated,
                                                           size_t numElements) {
#ifdef ACL_DEVICE
  assert( timeDerivatives != nullptr );
  assert( timeEvaluated != nullptr );
  assert( time >= expansionPoint );
  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");
  static_assert(kernel::gpu_derivativeTaylorExpansion::TmpMaxMemRequiredInBytes == 0);

  kernel::gpu_derivativeTaylorExpansion intKrnl;
  intKrnl.numElements = numElements;
  intKrnl.I = timeEvaluated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = const_cast<const RealT **>(timeDerivatives);
    intKrnl.extraOffset_dQ(i) = m_derivativesOffsets[i];
  }

  // iterate over time derivatives
  const RealT deltaT = time - expansionPoint;
  intKrnl.power = 1.0;
  for(int derivative = 0; derivative < ConvergenceOrder; ++derivative) {
    intKrnl.streamPtr = device.api->getDefaultStream();
    intKrnl.execute(derivative);
    intKrnl.power *= deltaT / static_cast<RealT>(derivative + 1);
  }
#else
  assert(false && "no implementation provided");
#endif
}

void computeDerivativeTaylorExpansion(RealT time,
                                                     RealT expansionPoint,
                                                     RealT const*  timeDerivatives,
                                                     RealT timeEvaluated[tensor::Q::size()],
                                                     unsigned order) {
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
  for(unsigned derivative = order; derivative < ConvergenceOrder; ++derivative) {
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

unsigned int* getDerivativesOffsets() {
  return m_derivativesOffsets;
}


};

}

#endif

