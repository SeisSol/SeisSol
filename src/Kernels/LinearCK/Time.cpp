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

#include "Kernels/LinearCK/TimeBase.h"

#include "GravitationalFreeSurfaceBC.h"
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
#include "Kernels/MemoryOps.h"

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
    krnl.power(der) = krnl.power(der - 1) * timeStepWidth / static_cast<real>(der + 1);
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

void Time::evaluate(const real* coeffs,
                    const real* timeDerivatives,
                    real timeEvaluated[tensor::Q::size()]) {
  /*
   * assert alignments.
   */
  assert((reinterpret_cast<uintptr_t>(timeDerivatives)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(timeEvaluated)) % Alignment == 0);

  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");

  kernel::derivativeTaylorExpansion krnl;
  krnl.I = timeEvaluated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQ(i) = timeDerivatives + yateto::computeFamilySize<tensor::dQ>(1, i);
    krnl.power(i) = coeffs[i];
  }
  krnl.execute();
}

void Time::evaluateBatched(const real* coeffs,
                           const real** timeDerivatives,
                           real** timeIntegratedDofs,
                           unsigned numElements,
                           seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  assert(timeDerivatives != nullptr);
  assert(timeIntegratedDofs != nullptr);
  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");
  static_assert(kernel::gpu_derivativeTaylorExpansion::TmpMaxMemRequiredInBytes == 0);

#ifndef DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS
  kernel::gpu_derivativeTaylorExpansion krnl;
  krnl.numElements = numElements;
  krnl.I = timeIntegratedDofs;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQ(i) = const_cast<const real**>(timeDerivatives);
    krnl.extraOffset_dQ(i) = yateto::computeFamilySize<tensor::dQ>(1, i);
    krnl.power(i) = coeffs[i];
  }
  krnl.streamPtr = runtime.stream();
  krnl.execute();
#else
  seissol::kernels::time::aux::taylorSum(numElements,
                                         timeIntegratedDofs,
                                         const_cast<const real**>(timeDerivatives),
                                         coeffs,
                                         runtime.stream());
#endif
#else
  logError() << "No GPU implementation provided";
#endif
}

void Time::flopsEvaluate(long long& nonZeroFlops, long long& hardwareFlops) {
  nonZeroFlops = kernel::derivativeTaylorExpansion::NonZeroFlops;
  hardwareFlops = kernel::derivativeTaylorExpansion::HardwareFlops;
}

void Time::setGlobalData(const CompoundGlobalData& global) {}

} // namespace seissol::kernels::solver::linearck
