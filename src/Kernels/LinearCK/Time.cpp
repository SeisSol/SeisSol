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

#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "GravitationalFreeSurfaceBC.h"
#include <Alignment.h>
#include <DataTypes/ConditionalTable.h>
#include <Initializer/BasicTypedefs.h>
#include <Kernels/Interface.h>
#include <Kernels/LinearCK/Solver.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Layer.h>
#include <Parallel/Runtime/Stream.h>
#include <algorithm>
#include <cstdint>

#include "Kernels/Common.h"
#include "Kernels/MemoryOps.h"

#include <cassert>
#include <cstring>
#include <stdint.h>

#include <yateto.h>
#include <yateto/InitTools.h>

#include "utils/logger.h"

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#endif

#include "Memory/GlobalData.h"

GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)

namespace seissol::kernels::solver::linearck {

template <typename Cfg>
void Spacetime<Cfg>::setGlobalData(const GlobalData& global) {
  m_krnlPrototype.kDivMT = global.get<Cfg>().stiffnessMatricesTransposed;
  projectDerivativeToNodalBoundaryRotated.V3mTo2nFace = global.get<Cfg>().v3mTo2nFace;

#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);

  deviceKrnlPrototype.kDivMT = global.get<Cfg, Executor::Device>().stiffnessMatricesTransposed;
  deviceDerivativeToNodalBoundaryRotated.V3mTo2nFace =
      global.get<Cfg, Executor::Device>().v3mTo2nFace;
#endif
}

template <typename Cfg>
void Spacetime<Cfg>::computeAder(const real* coeffs,
                                 double timeStepWidth,
                                 LTS::Ref<Cfg>& data,
                                 LocalTmp<Cfg>& tmp,
                                 real timeIntegrated[tensor::I<Cfg>::size()],
                                 real* timeDerivatives,
                                 bool updateDisplacement) {

  assert(reinterpret_cast<uintptr_t>(data.template get<LTS::Dofs>()) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(timeIntegrated) % Alignment == 0);
  assert(timeDerivatives == nullptr ||
         reinterpret_cast<uintptr_t>(timeDerivatives) % Alignment == 0);

  // Only a small fraction of cells has the gravitational free surface boundary condition
  updateDisplacement &=
      std::any_of(std::begin(data.template get<LTS::CellInformation>().faceTypes),
                  std::end(data.template get<LTS::CellInformation>().faceTypes),
                  [](const FaceType f) { return f == FaceType::FreeSurfaceGravity; });

  alignas(PagesizeStack) real temporaryBuffer[Solver::DerivativesSize<Cfg>];
  auto* derivativesBuffer = (timeDerivatives != nullptr) ? timeDerivatives : temporaryBuffer;

  kernel::derivative<Cfg> krnl = m_krnlPrototype;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
    krnl.star(i) = data.template get<LTS::LocalIntegration>().starMatrices[i];
  }

  // Optional source term
  set_ET(krnl, get_ptr_sourceMatrix(data.template get<LTS::LocalIntegration>().specific));

  krnl.dQ(0) = const_cast<real*>(data.template get<LTS::Dofs>());
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
    krnl.dQ(i) = derivativesBuffer + yateto::computeFamilySize<tensor::dQ<Cfg>>(1, i);
  }

  krnl.I = timeIntegrated;
  // powers in the taylor-series expansion
  for (std::size_t der = 0; der < Cfg::ConvergenceOrder; ++der) {
    krnl.power(der) = coeffs[der];
  }

  if (updateDisplacement) {
    // First derivative if needed later in kernel
    std::copy_n(data.template get<LTS::Dofs>(), tensor::dQ<Cfg>::size(0), derivativesBuffer);
  } else if (timeDerivatives != nullptr) {
    // First derivative is not needed here but later
    // Hence stream it out
    streamstore(tensor::dQ<Cfg>::size(0), data.template get<LTS::Dofs>(), derivativesBuffer);
  }

  krnl.execute();

  // Do not compute it like this if at interface
  // Compute integrated displacement over time step if needed.
  if (updateDisplacement) {
    auto& bc = tmp.gravitationalFreeSurfaceBc;
    for (std::size_t face = 0; face < 4; ++face) {
      if (data.template get<LTS::FaceDisplacements>()[face] != nullptr &&
          data.template get<LTS::CellInformation>().faceTypes[face] ==
              FaceType::FreeSurfaceGravity) {
        bc.evaluate(face,
                    projectDerivativeToNodalBoundaryRotated,
                    data.template get<LTS::BoundaryMapping>()[face],
                    data.template get<LTS::FaceDisplacements>()[face],
                    tmp.nodalAvgDisplacements[face].data(),
                    *this,
                    derivativesBuffer,
                    timeStepWidth,
                    data.template get<LTS::Material>(),
                    data.template get<LTS::CellInformation>().faceTypes[face]);
      }
    }
  }
}

template <typename Cfg>
void Spacetime<Cfg>::computeBatchedAder(const real* coeffs,
                                        double timeStepWidth,
                                        LocalTmp<Cfg>& tmp,
                                        ConditionalPointersToRealsTable& dataTable,
                                        ConditionalMaterialTable& materialTable,
                                        bool updateDisplacement,
                                        seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  kernel::gpu_derivative<Cfg> derivativesKrnl = deviceKrnlPrototype;

  ConditionalKey timeVolumeKernelKey(KernelNames::Time || KernelNames::Volume);
  if (dataTable.find(timeVolumeKernelKey) != dataTable.end()) {
    auto& entry = dataTable[timeVolumeKernelKey];

    const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
    derivativesKrnl.numElements = numElements;
    derivativesKrnl.I = (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr();

    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
      derivativesKrnl.star(i) = const_cast<const real**>(
          (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
      derivativesKrnl.extraOffset_star(i) =
          SEISSOL_ARRAY_OFFSET(LocalIntegrationData, starMatrices, i);
    }

    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
      derivativesKrnl.dQ(i) = (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr();
      derivativesKrnl.extraOffset_dQ(i) = yateto::computeFamilySize<tensor::dQ<Cfg>>(1, i);
    }

    // stream dofs to the zero derivative
    device.algorithms.streamBatchedData(
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr()),
        (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr(),
        tensor::Q<Cfg>::Size,
        derivativesKrnl.numElements,
        runtime.stream());

    const auto maxTmpMem = yateto::getMaxTmpMemRequired(derivativesKrnl);
    auto tmpMem = runtime.memoryHandle<real>((maxTmpMem * numElements) / sizeof(real));

    for (std::size_t der = 0; der < Cfg::ConvergenceOrder; ++der) {
      derivativesKrnl.power(der) = coeffs[der];
    }
    derivativesKrnl.linearAllocator.initialize(tmpMem.get());
    derivativesKrnl.streamPtr = runtime.stream();
    derivativesKrnl.execute();
  }

  if (updateDisplacement) {
    auto& bc = tmp.gravitationalFreeSurfaceBc;
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
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

template <typename Cfg>
void Spacetime<Cfg>::flopsAder(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) {
  nonZeroFlops = kernel::derivative<Cfg>::NonZeroFlops;
  hardwareFlops = kernel::derivative<Cfg>::HardwareFlops;
}

template <typename Cfg>
std::uint64_t Spacetime<Cfg>::bytesAder() {
  std::uint64_t reals = 0;

  // DOFs load, tDOFs load, tDOFs write
  reals += tensor::Q<Cfg>::size() + 2 * tensor::I<Cfg>::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star<Cfg>>();

  /// \todo incorporate derivatives

  return reals * sizeof(real);
}

template <typename Cfg>
void Time<Cfg>::evaluate(const real* coeffs,
                         const real* timeDerivatives,
                         real timeEvaluated[tensor::Q<Cfg>::size()]) {
  /*
   * assert alignments.
   */
  assert((reinterpret_cast<uintptr_t>(timeDerivatives)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(timeEvaluated)) % Alignment == 0);

  static_assert(tensor::I<Cfg>::size() == tensor::Q<Cfg>::size(),
                "Sizes of tensors I and Q must match");

  kernel::derivativeTaylorExpansion<Cfg> krnl;
  krnl.I = timeEvaluated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
    krnl.dQ(i) = timeDerivatives + yateto::computeFamilySize<tensor::dQ<Cfg>>(1, i);
    krnl.power(i) = coeffs[i];
  }
  krnl.execute();
}

template <typename Cfg>
void Time<Cfg>::evaluateBatched(const real* coeffs,
                                const real** timeDerivatives,
                                real** timeIntegratedDofs,
                                std::size_t numElements,
                                seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  assert(timeDerivatives != nullptr);
  assert(timeIntegratedDofs != nullptr);
  static_assert(tensor::I<Cfg>::size() == tensor::Q<Cfg>::size(),
                "Sizes of tensors I and Q must match");
  static_assert(kernel::gpu_derivativeTaylorExpansion<Cfg>::TmpMaxMemRequiredInBytes == 0);

#ifndef DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS
  kernel::gpu_derivativeTaylorExpansion<Cfg> krnl;
  krnl.numElements = numElements;
  krnl.I = timeIntegratedDofs;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
    krnl.dQ(i) = const_cast<const real**>(timeDerivatives);
    krnl.extraOffset_dQ(i) = yateto::computeFamilySize<tensor::dQ<Cfg>>(1, i);
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

template <typename Cfg>
void Time<Cfg>::flopsEvaluate(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) {
  nonZeroFlops = kernel::derivativeTaylorExpansion<Cfg>::NonZeroFlops;
  hardwareFlops = kernel::derivativeTaylorExpansion<Cfg>::HardwareFlops;
}

template <typename Cfg>
void Time<Cfg>::setGlobalData(const GlobalData& global) {}

#define SEISSOL_CONFIGITER(cfg) template class Spacetime<cfg>;
#include "ConfigIncludeLinearCK.h"

#define SEISSOL_CONFIGITER(cfg) template class Time<cfg>;
#include "ConfigIncludeLinearCK.h"

} // namespace seissol::kernels::solver::linearck
