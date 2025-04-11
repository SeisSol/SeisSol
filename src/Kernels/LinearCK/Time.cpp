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

#include "Kernels/Time.h"
#include "GravitationalFreeSurfaceBC.h"
#include "TimeBase.h"
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

GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)

namespace seissol::kernels::solver::linearck {
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

  alignas(PagesizeStack) real temporaryBuffer[yateto::computeFamilySize<tensor::dQ>()];
  auto* derivativesBuffer = (timeDerivatives != nullptr) ? timeDerivatives : temporaryBuffer;

  kernel::derivative krnl = m_krnlPrototype;
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
  krnl.power(0) = timeStepWidth;
  for (std::size_t der = 1; der < ConvergenceOrder; ++der) {
    krnl.power(der) = krnl.power(der - 1) * timeStepWidth / real(der + 1);
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

    unsigned starOffset = 0;
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      derivativesKrnl.star(i) =
          const_cast<const real**>((entry.get(inner_keys::Wp::Id::Star))->getDeviceDataPtr());
      derivativesKrnl.extraOffset_star(i) = starOffset;
      starOffset += tensor::star::size(i);
    }

    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
      derivativesKrnl.dQ(i) = (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr();
      derivativesKrnl.extraOffset_dQ(i) = yateto::computeFamilySize<tensor::dQ>(1, i);
    }

    // stream dofs to the zero derivative
    device.algorithms.streamBatchedData(
        (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr(),
        (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr(),
        tensor::Q::Size,
        derivativesKrnl.numElements,
        runtime.stream());

    const auto maxTmpMem = yateto::getMaxTmpMemRequired(derivativesKrnl);
    real* tmpMem = reinterpret_cast<real*>(device.api->getStackMemory(maxTmpMem * numElements));

    derivativesKrnl.power(0) = timeStepWidth;
    for (std::size_t Der = 1; Der < ConvergenceOrder; ++Der) {
      derivativesKrnl.power(Der) = derivativesKrnl.power(Der - 1) * timeStepWidth / real(Der + 1);
    }
    derivativesKrnl.linearAllocator.initialize(tmpMem);
    derivativesKrnl.streamPtr = runtime.stream();
    derivativesKrnl.execute();
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
  real* tmpMem = reinterpret_cast<real*>(
      device.api->getStackMemory(intKrnl.TmpMaxMemRequiredInBytes * numElements));

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
  intKrnl.linearAllocator.initialize(tmpMem);
  intKrnl.streamPtr = runtime.stream();
  intKrnl.execute();
  device.api->popStackMemory();
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
    intKrnl.power(derivative) = intKrnl.power(derivative - 1) * deltaT / real(derivative);
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

} // namespace seissol::kernels::solver::linearck
