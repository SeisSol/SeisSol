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

#include "Kernels/LinearCK/Time.h"

#include "Alignment.h"
#include "Common/Constants.h"
#include "Common/Marker.h"
#include "Config.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Interface.h"
#include "Kernels/LinearCK/Solver.h"
#include "Kernels/MemoryOps.h"
#include "Memory/Descriptor/LTS.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <stdint.h>
#include <utils/logger.h>
#include <yateto.h>
#include <yateto/InitTools.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#endif

GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(extraOffset_ET)
GENERATE_HAS_MEMBER(sourceMatrix)

namespace seissol::kernels::solver::linearck {
template <typename Cfg>
void Spacetime<Cfg>::setGlobalData(const CompoundGlobalData<Cfg>& global) {
  krnlPrototype_.bindGlobals(*global.onHost);
  fsgKernelPrototype_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);

  deviceKrnlPrototype_.bindGlobals(*global.onDevice);
  deviceFsgKernelPrototype_.bindGlobals(*global.onDevice);
#endif
}

template <typename Cfg>
void Spacetime<Cfg>::computeAder(const real* coeffs,
                                 double timeStepWidth,
                                 LTS::Ref<Cfg>& data,
                                 LocalTmp<Cfg>& tmp,
                                 real* timeIntegrated,
                                 real* timeDerivatives,
                                 bool updateDisplacement) {

  assert(reinterpret_cast<uintptr_t>(data.template get<LTS::Dofs>()) % Vectorsize == 0);
  assert(reinterpret_cast<uintptr_t>(timeIntegrated) % Vectorsize == 0);
  assert(timeDerivatives == nullptr ||
         reinterpret_cast<uintptr_t>(timeDerivatives) % Vectorsize == 0);

  // Only a small fraction of cells has the gravitational free surface boundary condition
  updateDisplacement &= [&]() {
    bool anyOfResult = false;
    for (std::size_t i = 0; i < Cell::NumFaces; ++i) {
      anyOfResult |=
          data.template get<LTS::CellInformation>().faceTypes[i] == FaceType::FreeSurfaceGravity;
    }
    return anyOfResult;
  }();

  alignas(PagesizeStack) real temporaryBuffer[Solver<Cfg>::DerivativesSize];
  auto* derivativesBuffer = (timeDerivatives != nullptr) ? timeDerivatives : temporaryBuffer;

  kernel::derivative<Cfg> krnl = krnlPrototype_;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
    krnl.star(i) = data.template get<LTS::LocalIntegration>().starMatrices[i];
  }

  // Optional source term
  set_ET(krnl, get_ptr_sourceMatrix(data.template get<LTS::LocalIntegration>().specific));

  krnl.dQ(0) = const_cast<real*>(data.template get<LTS::Dofs>());
  for (std::size_t i = 1; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
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
                    fsgKernelPrototype_,
                    data.template get<LTS::BoundaryMapping>()[face],
                    data.template get<LTS::FaceDisplacements>()[face],
                    tmp.nodalAvgDisplacements[face].data(),
                    derivativesBuffer,
                    coeffs,
                    timeStepWidth,
                    data.template get<LTS::Material>());
      }
    }
  }
}

template <typename Cfg>
void Spacetime<Cfg>::computeBatchedAder(
    SEISSOL_GPU_PARAM const real* coeffs,
    SEISSOL_GPU_PARAM double timeStepWidth,
    SEISSOL_GPU_PARAM LTS::Layer& layer,
    SEISSOL_GPU_PARAM LocalTmp<Cfg>& tmp,
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& dataTable,
    SEISSOL_GPU_PARAM bool updateDisplacement,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  using namespace seissol::recording;
  kernel::gpu_derivative<Cfg> derivativesKrnl = deviceKrnlPrototype_;

  const ConditionalKey timeVolumeKernelKey(KernelNames::Time || KernelNames::Volume);
  if (dataTable.find(timeVolumeKernelKey) != dataTable.end()) {
    auto& entry = dataTable[timeVolumeKernelKey];

    const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
    derivativesKrnl.numElements = numElements;
    derivativesKrnl.I = (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr();

    const auto** localIntegrationPtrs = const_cast<const real**>(
        (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());

    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData<Cfg>, starMatrices);
    for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
      derivativesKrnl.star(i) = localIntegrationPtrs;
      derivativesKrnl.extraOffset_star(i) =
          SEISSOL_ARRAY_OFFSET(LocalIntegrationData<Cfg>, starMatrices, i);
    }

    constexpr auto SourceMatrixOffset =
        offsetof(LocalIntegrationData<Cfg>, specific) +
        get_offset_sourceMatrix<decltype(LocalIntegrationData<Cfg>::specific)>();
    static_assert(SourceMatrixOffset % sizeof(real) == 0,
                  "SourceMatrixOffset is not dividable by the real size.");

    set_ET(derivativesKrnl, localIntegrationPtrs);
    set_extraOffset_ET(derivativesKrnl, SourceMatrixOffset / sizeof(real));

    for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
      derivativesKrnl.dQ(i) = (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr();
      derivativesKrnl.extraOffset_dQ(i) = yateto::computeFamilySize<tensor::dQ<Cfg>>(1, i);
    }

    derivativesKrnl.Q =
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr());

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
      bc.evaluateOnDevice(
          face, deviceFsgKernelPrototype_, dataTable, timeStepWidth, device_, runtime);
    }
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

template <typename Cfg>
PerformanceEstimate Spacetime<Cfg>::metrics() const {
  auto estimate = PerformanceEstimate::fromKernel<kernel::derivative<Cfg>>();

  // legacy memory estimate
  std::uint64_t reals = 0;

  // DOFs load, tDOFs load, tDOFs write
  reals += tensor::Q<Cfg>::size() + 2 * tensor::I<Cfg>::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star<Cfg>>();

  /// \todo incorporate derivatives

  estimate.bytes = reals * sizeof(real);

  return estimate;
}

template <typename Cfg>
void Time<Cfg>::evaluate(const real* coeffs,
                         const real* timeDerivatives,
                         real timeEvaluated[tensor::Q<Cfg>::size()]) {
  /*
   * assert alignments.
   */
  assert((reinterpret_cast<uintptr_t>(timeDerivatives)) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(timeEvaluated)) % Vectorsize == 0);

  static_assert(tensor::I<Cfg>::size() == tensor::Q<Cfg>::size(),
                "Sizes of tensors I and Q must match");

  kernel::derivativeTaylorExpansion<Cfg> krnl;
  krnl.I = timeEvaluated;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
    krnl.dQ(i) = timeDerivatives + yateto::computeFamilySize<tensor::dQ<Cfg>>(1, i);
    krnl.power(i) = coeffs[i];
  }
  krnl.execute();
}

template <typename Cfg>
void Time<Cfg>::evaluateBatched(
    SEISSOL_GPU_PARAM const real* coeffs,
    SEISSOL_GPU_PARAM const real** timeDerivatives,
    SEISSOL_GPU_PARAM real** timeIntegratedDofs,
    SEISSOL_GPU_PARAM std::size_t numElements,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;

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
    krnl.dQ(i) = timeDerivatives;
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
PerformanceEstimate Time<Cfg>::metrics() const {
  return PerformanceEstimate::fromKernel<kernel::derivativeTaylorExpansion<Cfg>>();
}

template <typename Cfg>
void Time<Cfg>::setGlobalData(const CompoundGlobalData<Cfg>& global) {}

#define SEISSOL_INSTANTIATE(Cfg)                                                                   \
  template class Spacetime<Cfg>;                                                                   \
  template class Time<Cfg>;
SEISSOL_FOR_EACH_CONFIG_LINEARCK(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels::solver::linearck
