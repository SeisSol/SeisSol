#include "Kernels/TimeBase.h"
#include "Kernels/Time.h"

#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

#include <Kernels/common.hpp>
#include <Kernels/denseMatrixOps.hpp>

#include <cstring>
#include <cassert>
#include <stdint.h>
#include <omp.h>
#include <Eigen/Dense>

#include "Equations/poroelastic/Model/PoroelasticSetup.h"

#include <yateto.h>

GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)

seissol::kernels::TimeBase::TimeBase(){
  m_derivativesOffsets[0] = 0;
  for (int order = 0; order < CONVERGENCE_ORDER; ++order) {
    if (order > 0) {
      m_derivativesOffsets[order] = tensor::dQ::size(order-1) + m_derivativesOffsets[order-1];
    }
  }
}

void seissol::kernels::Time::setHostGlobalData(GlobalData const* global) {
  for (int n = 0; n < CONVERGENCE_ORDER; ++n) {
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

void seissol::kernels::Time::setGlobalData(const CompoundGlobalData& global) {
  setHostGlobalData(global.onHost);

#ifdef ACL_DEVICE
  logError() << "Poroelasticity does not work on GPUs.";
#endif
}

void seissol::kernels::Time::executeSTP( double                      i_timeStepWidth,
                                         LocalData&                  data,
                                         real                        o_timeIntegrated[tensor::I::size()],
                                         real*                       stp )

{
  alignas(PAGESIZE_STACK) real stpRhs[tensor::spaceTimePredictorRhs::size()];
  assert( ((uintptr_t)stp) % ALIGNMENT == 0);
  std::fill(std::begin(stpRhs), std::end(stpRhs), 0);
  std::fill(stp, stp + tensor::spaceTimePredictor::size(), 0);
  kernel::spaceTimePredictor krnl = m_krnlPrototype;
 
  //libxsmm can not generate GEMMs with alpha!=1. As a workaround we multiply the 
  //star matrices with dt before we execute the kernel.
  real A_values[init::star::size(0)];
  real B_values[init::star::size(1)];
  real C_values[init::star::size(2)];
  for (size_t i = 0; i < init::star::size(0); i++) {
    A_values[i] = i_timeStepWidth * data.localIntegration.starMatrices[0][i];
    B_values[i] = i_timeStepWidth * data.localIntegration.starMatrices[1][i];
    C_values[i] = i_timeStepWidth * data.localIntegration.starMatrices[2][i];
  }
  krnl.star(0) = A_values;
  krnl.star(1) = B_values;
  krnl.star(2) = C_values;

  krnl.Gk = data.localIntegration.specific.G[10] * i_timeStepWidth;
  krnl.Gl = data.localIntegration.specific.G[11] * i_timeStepWidth;
  krnl.Gm = data.localIntegration.specific.G[12] * i_timeStepWidth;

  krnl.Q = const_cast<real*>(data.dofs);
  krnl.I = o_timeIntegrated;
  krnl.timestep = i_timeStepWidth;
  krnl.spaceTimePredictor = stp;
  krnl.spaceTimePredictorRhs = stpRhs;

  //The matrix Zinv depends on the timestep
  //If the timestep is not as expected e.g. when approaching a sync point
  //we have to recalculate it
  if (i_timeStepWidth != data.localIntegration.specific.typicalTimeStepWidth) {
    auto sourceMatrix = init::ET::view::create(data.localIntegration.specific.sourceMatrix);
    real ZinvData[NUMBER_OF_QUANTITIES][CONVERGENCE_ORDER*CONVERGENCE_ORDER];
    model::zInvInitializerForLoop<0, NUMBER_OF_QUANTITIES, decltype(sourceMatrix)>(ZinvData, sourceMatrix, i_timeStepWidth);
    for (size_t i = 0; i < NUMBER_OF_QUANTITIES; i++) {
      krnl.Zinv(i) = ZinvData[i];
    }
    // krnl.execute has to be run here: ZinvData is only allocated locally
    krnl.execute();
  } else {
    for (size_t i = 0; i < NUMBER_OF_QUANTITIES; i++) {
      krnl.Zinv(i) = data.localIntegration.specific.Zinv[i];
    }
    krnl.execute();
  }
}
                                          

void seissol::kernels::Time::computeAder( double i_timeStepWidth,
                                          LocalData& data,
                                          LocalTmp& tmp,
                                          real o_timeIntegrated[tensor::I::size()],
                                          real* o_timeDerivatives,
                                          bool updateDisplacement)
{
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)data.dofs)              % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeIntegrated )      % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeDerivatives)      % ALIGNMENT == 0 || o_timeDerivatives == NULL );

  alignas(ALIGNMENT) real temporaryBuffer[tensor::spaceTimePredictor::size()];
  real* stpBuffer = (o_timeDerivatives != nullptr) ? o_timeDerivatives : temporaryBuffer;
  executeSTP( i_timeStepWidth, data, o_timeIntegrated, stpBuffer );
}

void seissol::kernels::Time::evaluateAtTime(std::shared_ptr<seissol::basisFunction::SampledTimeBasisFunctions<real>> evaluatedTimeBasisFunctions,
                                            real const* timeDerivatives, real timeEvaluated[tensor::Q::size()]) {
  kernel::evaluateDOFSAtTimeSTP krnl;
  krnl.spaceTimePredictor = timeDerivatives;
  krnl.QAtTimeSTP = timeEvaluated;
  krnl.timeBasisFunctionsAtPoint = evaluatedTimeBasisFunctions->m_data.data();
  krnl.execute();
}

void flopsEvaluateAtTime(long long& nonZeroFlops, long long& hardwareFlops) {
  // reset flops
  nonZeroFlops = 0; hardwareFlops = 0;

  nonZeroFlops  += kernel::evaluateDOFSAtTimeSTP::NonZeroFlops;
  hardwareFlops += kernel::evaluateDOFSAtTimeSTP::HardwareFlops;
}

void seissol::kernels::Time::flopsAder( unsigned int        &o_nonZeroFlops,
                                        unsigned int        &o_hardwareFlops ) {
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops =0;

  o_nonZeroFlops = kernel::spaceTimePredictor::NonZeroFlops;
  o_hardwareFlops = kernel::spaceTimePredictor::HardwareFlops;
  //we multiply the star matrices with dt before we execute the kernel
  o_nonZeroFlops += 3*init::star::size(0);
  o_hardwareFlops += 3*init::star::size(0);
}

unsigned seissol::kernels::Time::bytesAder()
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

  return reals * sizeof(real);
}

void seissol::kernels::Time::computeIntegral( double                            i_expansionPoint,
                                              double                            i_integrationStart,
                                              double                            i_integrationEnd,
                                              const real*                       i_timeDerivatives,
                                              real                              o_timeIntegrated[tensor::I::size()])
{
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)i_timeDerivatives)  % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeIntegrated)   % ALIGNMENT == 0 );

  // assert that this is a forwared integration in time
  assert( i_integrationStart + (real) 1.E-10 > i_expansionPoint   );
  assert( i_integrationEnd                   > i_integrationStart );

  /*
   * compute time integral.
   */
  // compute lengths of integration intervals
  real l_deltaTLower = i_integrationStart - i_expansionPoint;
  real l_deltaTUpper = i_integrationEnd   - i_expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  real l_firstTerm  = (real) 1;
  real l_secondTerm = (real) 1;
  real l_factorial  = (real) 1;
  
  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = o_timeIntegrated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = i_timeDerivatives + m_derivativesOffsets[i];
  }
 
  // iterate over time derivatives
  for(int der = 0; der < CONVERGENCE_ORDER; ++der ) {
    l_firstTerm  *= l_deltaTUpper;
    l_secondTerm *= l_deltaTLower;
    l_factorial  *= (real)(der+1);

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

  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = timeEvaluated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = timeDerivatives + m_derivativesOffsets[i];
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
    nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(der);
    hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(der);
  }
}

void seissol::kernels::Time::computeBatchedIntegral(double i_expansionPoint,
                                                    double i_integrationStart,
                                                    double i_integrationEnd,
                                                    const real** i_timeDerivatives,
                                                    real ** o_timeIntegratedDofs,
                                                    unsigned numElements) {
#ifdef ACL_DEVICE
  // assert that this is a forwared integration in time
  assert( i_integrationStart + (real) 1.E-10 > i_expansionPoint   );
  assert( i_integrationEnd                   > i_integrationStart );

  /*
   * compute time integral.
   */
  // compute lengths of integration intervals
  real deltaTLower = i_integrationStart - i_expansionPoint;
  real deltaTUpper = i_integrationEnd - i_expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  real firstTerm  = static_cast<real>(1.0);
  real secondTerm = static_cast<real>(1.0);
  real factorial  = static_cast<real>(1.0);

  kernel::gpu_derivativeTaylorExpansion intKrnl;
  intKrnl.numElements = numElements;
  real* tmpMem = reinterpret_cast<real*>(device.api->getStackMemory(intKrnl.TmpMaxMemRequiredInBytes * numElements));

  intKrnl.I = o_timeIntegratedDofs;

  unsigned derivativesOffset = 0;
  for (size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = i_timeDerivatives;
    intKrnl.extraOffset_dQ(i) = derivativesOffset;
    derivativesOffset += tensor::dQ::size(i);
  }

  // iterate over time derivatives
  for(int der = 0; der < CONVERGENCE_ORDER; ++der) {
    firstTerm *= deltaTUpper;
    secondTerm *= deltaTLower;
    factorial *= static_cast<real>(der + 1);

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

void seissol::kernels::Time::computeBatchedTaylorExpansion(real time,
                                                           real expansionPoint,
                                                           real** timeDerivatives,
                                                           real** timeEvaluated,
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
    intKrnl.dQ(i) = const_cast<const real **>(timeDerivatives);
    intKrnl.extraOffset_dQ(i) = m_derivativesOffsets[i];
  }

  // iterate over time derivatives
  const real deltaT = time - expansionPoint;
  intKrnl.power = 1.0;
  for(int derivative = 0; derivative < CONVERGENCE_ORDER; ++derivative) {
    intKrnl.streamPtr = device.api->getDefaultStream();
    intKrnl.execute(derivative);
    intKrnl.power *= deltaT / static_cast<real>(derivative + 1);
  }
#else
  assert(false && "no implementation provided");
#endif
}

void seissol::kernels::Time::computeBatchedAder(double i_timeStepWidth,
                                                LocalTmp& tmp,
                                                ConditionalPointersToRealsTable &dataTable,
                                                ConditionalMaterialTable &materialTable,
                                                bool updateDisplacement) {
#ifdef ACL_DEVICE
  alignas(PAGESIZE_STACK) real stpRhs[tensor::spaceTimePredictorRhs::size()];
  assert( ((uintptr_t)stp) % ALIGNMENT == 0);
  std::fill(std::begin(stpRhs), std::end(stpRhs), 0);
  std::fill(stp, stp + tensor::spaceTimePredictor::size(), 0);
  kernel::gpu_spaceTimePredictor krnl;

  ConditionalKey timeVolumeKernelKey(KernelNames::Time || KernelNames::Volume);
  if(dataTable.find(timeVolumeKernelKey) != dataTable.end()) {
    auto &entry = dataTable[timeVolumeKernelKey];

    const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
    krnl.numElements = numElements;

    krnl.I = (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr();
    krnl.Q = (entry.get(inner_keys::Wp::Id::Qdofs))->getDeviceDataPtr();
    krnl.timestep = i_timeStepWidth;
    krnl.spaceTimePredictor = (entry.get(inner_keys::Wp::Id::Stp))->getDeviceDataPtr();
    krnl.spaceTimePredictorRhs = (entry.get(inner_keys::Wp::Id::StpRhs))->getDeviceDataPtr();

    std::size_t starOffset = 0;
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      krnl.star(i) = const_cast<const real **>((entry.get(inner_keys::Wp::Id::Star))->getDeviceDataPtr());
      krnl.extraOffset_star(i) = starOffset;
      starOffset += tensor::star::size(i);
    }

    krnl.Gk = data.localIntegration.specific.G[10] * i_timeStepWidth;
    krnl.Gl = data.localIntegration.specific.G[11] * i_timeStepWidth;
    krnl.Gm = data.localIntegration.specific.G[12] * i_timeStepWidth;

    if (i_timeStepWidth != data.localIntegration.specific.typicalTimeStepWidth) {
      assert(false && "NYI");
    }
    else {
      std::size_t zinvOffset = 0;
      for (size_t i = 0; i < yateto::numFamilyMembers<tensor::Zinv>(); i++) {
        krnl.Zinv(i) = const_cast<const real **>((entry.get(inner_keys::Wp::Id::Zinv))->getDeviceDataPtr());
        krnl.extraOffset_Zinv(i) = zinvOffset;
        zinvOffset += tensor::Zinv::size(i);
      }
      krnl.execute();
    }
  }
#else
  assert(false && "no implementation provided");
#endif
}
