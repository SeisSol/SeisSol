// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#include "Kernels/NonLinearCK/TimeBase.h"

#include "Kernels/LinearCK/GravitationalFreeSurfaceBC.h"
#include <Alignment.h>
#include <Common/Constants.h>
#include <DataTypes/ConditionalTable.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Interface.h>
#include <Kernels/Precision.h>
#include <Numerical/BasisFunction.h>
#include <Parallel/Runtime/Stream.h>
#include <algorithm>
#include <generated_code/kernel.h>
#include <generated_code/tensor.h>
#include <iterator>

#include "Kernels/Common.h"
#include "Kernels/DenseMatrixOps.h"

#include <cassert>
#include <cstring>
#include <memory>
#include <stdint.h>

#include <yateto.h>
#include <yateto/InitTools.h>

#include "utils/logger.h"

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#endif

GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)

namespace seissol::kernels::solver::nonlinearck {
void Spacetime::setGlobalData(const CompoundGlobalData& global) {
  m_krnlPrototype.kDivMT = global.onHost->stiffnessMatricesTransposed;
  projectDerivativeToNodalBoundaryRotated.V3mTo2nFace = global.onHost->V3mTo2nFace;

#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);
  const auto deviceAlignment = device.api->getGlobMemAlignment();
  deviceKrnlPrototype.kDivMT = global.onDevice->stiffnessMatricesTransposed;
  deviceDerivativeToNodalBoundaryRotated.V3mTo2nFace = global.onDevice->V3mTo2nFace;
#endif
}

void Spacetime::computeAder(double timeStepWidth,
                       LocalData& data,
                       LocalTmp& tmp,
                       real timeIntegrated[tensor::I::size()],
                       real* timeDerivatives,
                       bool updateDisplacement) {

  assert(reinterpret_cast<uintptr_t>(data.dofs()) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(timeIntegrated) % Alignment == 0);
  assert(timeDerivatives == nullptr ||
         reinterpret_cast<uintptr_t>(timeDerivatives) % Alignment == 0);

  // Only a small fraction of cells has the gravitational free surface boundary condition
  updateDisplacement &=
      std::any_of(std::begin(data.cellInformation().faceTypes),
                  std::end(data.cellInformation().faceTypes),
                  [](const FaceType f) { return f == FaceType::FreeSurfaceGravity; });

// // In this computeAder() function, 'I' and 'derivatives in linear case is computed;
// // In nonlinear case, 'derivatives' and 'F' will be computed here.
// // Step 1: Compute each order of derivatives
//   alignas(PagesizeStack) real temporaryBuffer[yateto::computeFamilySize<tensor::dQ>()];
//   auto* derivativesBuffer = (timeDerivatives != nullptr) ? timeDerivatives : temporaryBuffer;

//   // Step 1.1: Convert the Modal solution data.dofs(), Q, to Nodal space
//   kernel::damageConvertToNodal d_converToKrnl;
//   alignas(PagesizeStack) real solNData[tensor::QNodal::size()];
//   d_converToKrnl.v = init::v::Values;
//   d_converToKrnl.QNodal = solNData;
//   d_converToKrnl.Q = data.dofs();
//   d_converToKrnl.execute();

//   // Step 1.2: Compute rhs of damage evolution
//   alignas(PagesizeStack) real fNodalData[tensor::FNodal::size()] = {0};
//   real* exxNodal = ( solNData + 0*tensor::Q::Shape[0] );
//   real* eyyNodal = (solNData + 1*tensor::Q::Shape[0]);
//   real* ezzNodal = (solNData + 2*tensor::Q::Shape[0]);
//   real* exyNodal = (solNData + 3*tensor::Q::Shape[0]);
//   real* eyzNodal = (solNData + 4*tensor::Q::Shape[0]);
//   real* ezxNodal = (solNData + 5*tensor::Q::Shape[0]);
//   real* alphaNodal = (solNData + 9*tensor::Q::Shape[0]);
//   // real* breakNodal = (solNData + 10*tensor::Q::Shape[0]);

//   real alpha_ave = 0.0;
//   // real break_ave = 0.0;
//   real w_ave = 1.0/tensor::Q::Shape[0];
//   for (unsigned int q = 0; q<tensor::Q::Shape[0]; ++q){
//     // break_ave += breakNodal[q] * w_ave;
//     alpha_ave += alphaNodal[q] * w_ave;
//   }

//   real const damage_para1 = data.material().local.Cd; // 1.2e-4*2;
//   real const damage_para2 = data.material().local.gammaR;

//   for (unsigned int q = 0; q<tensor::Q::Shape[0]; ++q){
//     real EspI = (exxNodal[q]+data.material().local.epsInit_xx) + 
//       (eyyNodal[q]+data.material().local.epsInit_yy) + (ezzNodal[q]+data.material().local.epsInit_zz);
//     real EspII = 
//       (exxNodal[q]+data.material().local.epsInit_xx)*(exxNodal[q]+data.material().local.epsInit_xx)
//     + (eyyNodal[q]+data.material().local.epsInit_yy)*(eyyNodal[q]+data.material().local.epsInit_yy)
//     + (ezzNodal[q]+data.material().local.epsInit_zz)*(ezzNodal[q]+data.material().local.epsInit_zz)
//     + 2*(exyNodal[q]+data.material().local.epsInit_xy)*(exyNodal[q]+data.material().local.epsInit_xy)
//     + 2*(eyzNodal[q]+data.material().local.epsInit_yz)*(eyzNodal[q]+data.material().local.epsInit_yz)
//     + 2*(ezxNodal[q]+data.material().local.epsInit_xz)*(ezxNodal[q]+data.material().local.epsInit_xz);

//     real W_energy = 0.0*0.5*data.material().local.lambdaE*EspI*EspI
//         + data.material().local.muE*EspII;

//     if (W_energy - damage_para2*(alpha_ave/(1-alpha_ave))*(alpha_ave/(1-alpha_ave)) > 0) {
//       if (alpha_ave < 0.8){
//         fNodalData[9*tensor::Q::Shape[0] + q] =
//           1.0/(damage_para1*damage_para2)
//                 *(W_energy - damage_para2*(alpha_ave/(1-alpha_ave))*(alpha_ave/(1-alpha_ave)));
//       }
//       else{
//         fNodalData[9*tensor::Q::Shape[0] + q] = 0.0;
//       }
//     } else if (alpha_ave > 8e-1) {
//       fNodalData[9*tensor::Q::Shape[0] + q] =
//         1.0/(damage_para1*damage_para2)
//                 *(W_energy - damage_para2*(alpha_ave/(1-alpha_ave))*(alpha_ave/(1-alpha_ave)));
//     }
//     else {
//       fNodalData[9*tensor::Q::Shape[0] + q] = 0;
//     }

//   }

// // Step 2: Convert from Modal to Nodal for each temporal quadrature point;
// // Meanwhile, compute the nonlinear nodal Rusanov fluxes 
// // (TimeCluster.h, in previous version)
// // For neighbor cell: 0.5(Fpd * nd - C * up)
// // NOTE: need to include face relation in the final integration, but maybe 
// // not necessarily here - can be integrated in Neighbor.cpp

// // Step 3: Do time integration for the Rusanov flux

// // Step 4: Convert the integrated Rusanov flux from Nodal to Modal space

  alignas(PagesizeStack) real temporaryBuffer[yateto::computeFamilySize<tensor::dQ>()];
  auto* derivativesBuffer = (timeDerivatives != nullptr) ? timeDerivatives : temporaryBuffer;

  kernel::derivativeDamage krnl = m_krnlPrototype;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    krnl.star(i) = data.localIntegration().starMatrices[i];
  }

  // Optional source term
  set_ET(krnl, get_ptr_sourceMatrix(data.localIntegration().specific));

  krnl.dQ(0) = const_cast<real*>(data.dofs());
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQ(i) = derivativesBuffer + yateto::computeFamilySize<tensor::dQ>(1, i);
  }

  krnl.I = timeIntegrated;
  // powers in the taylor-series expansion
  krnl.nlPower(0) = timeStepWidth;
  for (std::size_t der = 1; der < ConvergenceOrder; ++der) {
    krnl.nlPower(der) = krnl.nlPower(der - 1) * timeStepWidth / real(der + 1);
  }

  if (updateDisplacement) {
    // First derivative if needed later in kernel
    std::copy_n(data.dofs(), tensor::dQ::size(0), derivativesBuffer);
  } else if (timeDerivatives != nullptr) {
    // First derivative is not needed here but later
    // Hence stream it out
    streamstore(tensor::dQ::size(0), data.dofs(), derivativesBuffer);
  }

  krnl.execute();

  // Do not compute it like this if at interface
  // Compute integrated displacement over time step if needed.
  if (updateDisplacement) {
    auto& bc = tmp.gravitationalFreeSurfaceBc;
    for (unsigned face = 0; face < 4; ++face) {
      if (data.faceDisplacements()[face] != nullptr &&
          data.cellInformation().faceTypes[face] == FaceType::FreeSurfaceGravity) {
        bc.evaluate(face,
                    projectDerivativeToNodalBoundaryRotated,
                    data.boundaryMapping()[face],
                    data.faceDisplacements()[face],
                    tmp.nodalAvgDisplacements[face].data(),
                    *this,
                    derivativesBuffer,
                    timeStepWidth,
                    data.material(),
                    data.cellInformation().faceTypes[face]);
      }
    }
  }
}

void Spacetime::computeBatchedAder(double timeStepWidth,
                                  LocalTmp& tmp,
                                  ConditionalPointersToRealsTable& dataTable,
                                  ConditionalMaterialTable& materialTable,
                                  bool updateDisplacement,
                                  seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  kernel::gpu_derivative derivativesKrnl = deviceKrnlPrototype;

  ConditionalKey timeVolumeKernelKey(KernelNames::Time || KernelNames::Volume);
  if (dataTable.find(timeVolumeKernelKey) != dataTable.end()) {
    auto& entry = dataTable[timeVolumeKernelKey];

    const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
    derivativesKrnl.numElements = numElements;
    derivativesKrnl.I = (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr();

    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      derivativesKrnl.star(i) = const_cast<const real**>(
          (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
      derivativesKrnl.extraOffset_star(i) =
          SEISSOL_ARRAY_OFFSET(LocalIntegrationData, starMatrices, i);
    }

    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
      derivativesKrnl.dQ(i) = (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr();
      derivativesKrnl.extraOffset_dQ(i) = yateto::computeFamilySize<tensor::dQ>(1, i);
    }

    // stream dofs to the zero derivative
    device.algorithms.streamBatchedData(
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr()),
        (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr(),
        tensor::Q::Size,
        derivativesKrnl.numElements,
        runtime.stream());

    const auto maxTmpMem = yateto::getMaxTmpMemRequired(derivativesKrnl);
    auto tmpMem = runtime.memoryHandle<real>((maxTmpMem * numElements) / sizeof(real));

    derivativesKrnl.power(0) = timeStepWidth;
    for (std::size_t Der = 1; Der < ConvergenceOrder; ++Der) {
      derivativesKrnl.power(Der) = derivativesKrnl.power(Der - 1) * timeStepWidth / real(Der + 1);
    }
    derivativesKrnl.linearAllocator.initialize(tmpMem.get());
    derivativesKrnl.streamPtr = runtime.stream();
    derivativesKrnl.execute();
  }

  if (updateDisplacement) {
    auto& bc = tmp.gravitationalFreeSurfaceBc;
    for (unsigned face = 0; face < 4; ++face) {
      bc.evaluateOnDevice(face,
                          deviceDerivativeToNodalBoundaryRotated,
                          *this,
                          dataTable,
                          materialTable,
                          timeStepWidth,
                          device,
                          runtime);
    }
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

void Spacetime::flopsAder(unsigned int& nonZeroFlops, unsigned int& hardwareFlops) {
  nonZeroFlops = kernel::derivative::NonZeroFlops;
  hardwareFlops = kernel::derivative::HardwareFlops;
}

unsigned Spacetime::bytesAder() {
  unsigned reals = 0;

  // DOFs load, tDOFs load, tDOFs write
  reals += tensor::Q::size() + 2 * tensor::I::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star>();

  /// \todo incorporate derivatives

  return reals * sizeof(real);
}

void Time::computeIntegral(double expansionPoint,
                           double integrationStart,
                           double integrationEnd,
                           const real* timeDerivatives,
                           real timeIntegrated[tensor::I::size()]) {
  /*
   * assert alignments.
   */
  assert((reinterpret_cast<uintptr_t>(timeDerivatives)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(timeIntegrated)) % Alignment == 0);

  // assert that this is a forwared integration in time
  assert(integrationStart + (real)1.E-10 > expansionPoint);
  assert(integrationEnd > integrationStart);

  /*
   * compute time integral.
   */
  // compute lengths of integration intervals
  const real deltaTLower = integrationStart - expansionPoint;
  const real deltaTUpper = integrationEnd - expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  real firstTerm = (real)1;
  real secondTerm = (real)1;
  real factorial = (real)1;

  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = timeIntegrated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = timeDerivatives + yateto::computeFamilySize<tensor::dQ>(1, i);
  }

  // iterate over time derivatives
  for (std::size_t der = 0; der < ConvergenceOrder; ++der) {
    firstTerm *= deltaTUpper;
    secondTerm *= deltaTLower;
    factorial *= (real)(der + 1);

    intKrnl.power(der) = firstTerm - secondTerm;
    intKrnl.power(der) /= factorial;
  }
  intKrnl.execute();
}

void Time::computeBatchedIntegral(double expansionPoint,
                                  double integrationStart,
                                  double integrationEnd,
                                  const real** timeDerivatives,
                                  real** timeIntegratedDofs,
                                  unsigned numElements,
                                  seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  // assert that this is a forwared integration in time
  assert(integrationStart + (real)1.E-10 > expansionPoint);
  assert(integrationEnd > integrationStart);

  /*
   * compute time integral.
   */
  // compute lengths of integration intervals
  const real deltaTLower = integrationStart - expansionPoint;
  const real deltaTUpper = integrationEnd - expansionPoint;

#ifndef DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS
  // compute lengths of integration intervals

  // initialization of scalars in the taylor series expansion (0th term)
  real firstTerm = static_cast<real>(1.0);
  real secondTerm = static_cast<real>(1.0);
  real factorial = static_cast<real>(1.0);

  kernel::gpu_derivativeTaylorExpansion intKrnl;
  intKrnl.numElements = numElements;
  auto tmpMem =
      runtime.memoryHandle<real>((intKrnl.TmpMaxMemRequiredInBytes * numElements) / sizeof(real));

  intKrnl.I = timeIntegratedDofs;

  for (size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = timeDerivatives;
    intKrnl.extraOffset_dQ(i) = yateto::computeFamilySize<tensor::dQ>(1, i);
  }

  // iterate over time derivatives
  for (std::size_t der = 0; der < ConvergenceOrder; ++der) {
    firstTerm *= deltaTUpper;
    secondTerm *= deltaTLower;
    factorial *= static_cast<real>(der + 1);

    intKrnl.power(der) = firstTerm - secondTerm;
    intKrnl.power(der) /= factorial;
  }
  intKrnl.linearAllocator.initialize(tmpMem.get());
  intKrnl.streamPtr = runtime.stream();
  intKrnl.execute();
#else
  seissol::kernels::time::aux::taylorSum(true,
                                         numElements,
                                         timeIntegratedDofs,
                                         timeDerivatives,
                                         deltaTLower,
                                         deltaTUpper,
                                         runtime.stream());
#endif
#else
  logError() << "No GPU implementation provided";
#endif
}

void Time::computeTaylorExpansion(real time,
                                  real expansionPoint,
                                  const real* timeDerivatives,
                                  real timeEvaluated[tensor::Q::size()]) {
  /*
   * assert alignments.
   */
  assert((reinterpret_cast<uintptr_t>(timeDerivatives)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(timeEvaluated)) % Alignment == 0);

  // assert that this is a forward evaluation in time
  assert(time >= expansionPoint);

  const real deltaT = time - expansionPoint;

  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");

  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = timeEvaluated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = timeDerivatives + yateto::computeFamilySize<tensor::dQ>(1, i);
  }
  intKrnl.power(0) = 1.0;

  // iterate over time derivatives
  for (std::size_t derivative = 1; derivative < ConvergenceOrder; ++derivative) {
    intKrnl.power(derivative) =
        intKrnl.power(derivative - 1) * deltaT / static_cast<real>(derivative);
  }

  intKrnl.execute();
}

void Time::computeBatchedTaylorExpansion(real time,
                                         real expansionPoint,
                                         real** timeDerivatives,
                                         real** timeEvaluated,
                                         size_t numElements,
                                         seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  assert(timeDerivatives != nullptr);
  assert(timeEvaluated != nullptr);
  assert(time >= expansionPoint);
  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");
  static_assert(kernel::gpu_derivativeTaylorExpansion::TmpMaxMemRequiredInBytes == 0);

  const real deltaT = time - expansionPoint;

#ifndef DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS
  kernel::gpu_derivativeTaylorExpansion intKrnl;
  intKrnl.numElements = numElements;
  intKrnl.I = timeEvaluated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = const_cast<const real**>(timeDerivatives);
    intKrnl.extraOffset_dQ(i) = yateto::computeFamilySize<tensor::dQ>(1, i);
  }

  // iterate over time derivatives
  intKrnl.power(0) = 1.0;
  for (std::size_t derivative = 1; derivative < ConvergenceOrder; ++derivative) {
    intKrnl.power(derivative) =
        intKrnl.power(derivative - 1) * deltaT / static_cast<real>(derivative);
  }

  intKrnl.streamPtr = runtime.stream();
  intKrnl.execute();
#else
  seissol::kernels::time::aux::taylorSum(false,
                                         numElements,
                                         timeEvaluated,
                                         const_cast<const real**>(timeDerivatives),
                                         0,
                                         deltaT,
                                         runtime.stream());
#endif
#else
  logError() << "No GPU implementation provided";
#endif
}

void Time::flopsTaylorExpansion(long long& nonZeroFlops, long long& hardwareFlops) {
  nonZeroFlops = kernel::derivativeTaylorExpansion::NonZeroFlops;
  hardwareFlops = kernel::derivativeTaylorExpansion::HardwareFlops;
}

void Time::evaluateAtTime(std::shared_ptr<seissol::basisFunction::SampledTimeBasisFunctions<real>>
                              evaluatedTimeBasisFunctions,
                          const real* timeDerivatives,
                          real timeEvaluated[tensor::Q::size()]) {
#ifdef USE_STP
  kernel::evaluateDOFSAtTimeSTP krnl;
  krnl.spaceTimePredictor = timeDerivatives;
  krnl.QAtTimeSTP = timeEvaluated;
  krnl.timeBasisFunctionsAtPoint = evaluatedTimeBasisFunctions->m_data.data();
  krnl.execute();
#endif
}

void Time::flopsEvaluateAtTime(long long& nonZeroFlops, long long& hardwareFlops) {
#ifdef USE_STP
  // reset flops
  nonZeroFlops = 0;
  hardwareFlops = 0;

  nonZeroFlops += kernel::evaluateDOFSAtTimeSTP::NonZeroFlops;
  hardwareFlops += kernel::evaluateDOFSAtTimeSTP::HardwareFlops;
#endif
}

void Time::setGlobalData(const CompoundGlobalData& global) {}

} // namespace seissol::kernels::solver::nonlinearck
