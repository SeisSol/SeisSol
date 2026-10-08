// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "Time.h"

#include "Alignment.h"
#include "Common/Constants.h"
#include "Common/Marker.h"
#include "Config.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Interface.h"
#include "Kernels/LinearCKAnelastic/Solver.h"
#include "Kernels/MemoryOps.h"
#include "Memory/Descriptor/LTS.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdint.h>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#endif

#ifndef ACL_DEVICE
#include <utils/logger.h>
#endif

namespace seissol::kernels::solver::linearckanelastic {

template <typename Cfg>
void Time<Cfg>::setGlobalData(const CompoundGlobalData<Cfg>& global) {}

template <typename Cfg>
void Spacetime<Cfg>::setGlobalData(const CompoundGlobalData<Cfg>& global) {
  krnlPrototype_.bindGlobals(*global.onHost);
  fsgKernelPrototype_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceKrnlPrototype_.bindGlobals(*global.onDevice);
#endif
}

template <typename Cfg>
void Spacetime<Cfg>::computeAder(const real* coeffs,
                                 double timeStepWidth,
                                 LTS::Ref<Cfg>& data,
                                 LocalTmp<Cfg>& tmp,
                                 real* timeIntegrated,
                                 real* timeDerivativesOrSTP,
                                 bool updateDisplacement) {
  /*
   * assert alignments.
   */
  assert((reinterpret_cast<uintptr_t>(data.template get<LTS::Dofs>())) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(timeIntegrated)) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(timeDerivativesOrSTP)) % Vectorsize == 0 ||
         timeDerivativesOrSTP == nullptr);

  /*
   * compute ADER scheme.
   */
  // temporary result
  // TODO(David): move these temporary buffers into the Yateto kernel, maybe
  alignas(PagesizeStack) real temporaryBuffer[2][tensor::dQ<Cfg>::size(0)];
  alignas(PagesizeStack) real temporaryBufferExt[2][tensor::dQext<Cfg>::size(1)];
  alignas(PagesizeStack) real temporaryBufferAne[2][tensor::dQane<Cfg>::size(0)];

  kernel::derivative<Cfg> krnl = krnlPrototype_;

  // Only a small fraction of cells has the gravitational free surface boundary condition
  updateDisplacement &= [&]() {
    bool anyOfResult = false;
    for (std::size_t i = 0; i < Cell::NumFaces; ++i) {
      anyOfResult |=
          data.template get<LTS::CellInformation>().faceTypes[i] == FaceType::FreeSurfaceGravity;
    }
    return anyOfResult;
  }();

  // the gravitational free surface boundary condition reads every derivative, so
  // they have to stay around for the whole timestep
  alignas(PagesizeStack) real derivativesScratch[Solver<Cfg>::DerivativesSize];
  real* derivativesBuffer = timeDerivativesOrSTP;
  if (derivativesBuffer == nullptr && updateDisplacement) {
    derivativesBuffer = derivativesScratch;
  }

  krnl.dQ(0) = const_cast<real*>(data.template get<LTS::Dofs>());
  if (derivativesBuffer != nullptr) {
    if (updateDisplacement) {
      std::copy_n(data.template get<LTS::Dofs>(), tensor::dQ<Cfg>::size(0), derivativesBuffer);
    } else {
      streamstore(tensor::dQ<Cfg>::size(0), data.template get<LTS::Dofs>(), derivativesBuffer);
    }
    for (std::size_t i = 1; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
      krnl.dQ(i) = derivativesBuffer + yateto::computeFamilySize<tensor::dQ<Cfg>>(1, i);
    }
  } else {
    for (std::size_t i = 1; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
      krnl.dQ(i) = temporaryBuffer[i % 2];
    }
  }

  krnl.dQane(0) = const_cast<real*>(data.template get<LTS::DofsAne>());
  for (std::size_t i = 1; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
    krnl.dQane(i) = temporaryBufferAne[i % 2];
    krnl.dQext(i) = temporaryBufferExt[i % 2];
  }

  krnl.I = timeIntegrated;
  krnl.Iane = tmp.timeIntegratedAne;

  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
    krnl.star(i) = data.template get<LTS::LocalIntegration>().starMatrices[i];
  }
  krnl.w = data.template get<LTS::LocalIntegration>().specific.w;
  krnl.W = data.template get<LTS::LocalIntegration>().specific.W;
  krnl.E = data.template get<LTS::LocalIntegration>().specific.E;

  // powers in the taylor-series expansion
  for (std::size_t der = 0; der < Cfg::ConvergenceOrder; ++der) {
    krnl.power(der) = coeffs[der];
  }

  krnl.execute();

  // Compute integrated displacement over time step if needed.
  if (updateDisplacement) {
    auto& bc = tmp.gravitationalFreeSurfaceBc;
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
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
PerformanceEstimate Spacetime<Cfg>::metrics() const {
  auto estimate = PerformanceEstimate::fromKernel<kernel::derivative<Cfg>>();

  // legacy memory estimate
  std::uint64_t reals = 0;

  // DOFs load, tDOFs load, tDOFs write
  reals += tensor::Q<Cfg>::size() + tensor::Qane<Cfg>::size() + 2 * tensor::I<Cfg>::size() +
           2 * tensor::Iane<Cfg>::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star<Cfg>>() + tensor::w<Cfg>::size() +
           tensor::W<Cfg>::size() + tensor::E<Cfg>::size();

  /// \todo incorporate derivatives

  estimate.bytes = reals * sizeof(real);

  return estimate;
}

template <typename Cfg>
void Time<Cfg>::evaluate(const real* coeffs,
                         const real* timeDerivativesOrSTP,
                         real timeEvaluated[tensor::Q<Cfg>::size()]) {
  /*
   * assert alignments.
   */
  assert((reinterpret_cast<uintptr_t>(timeDerivativesOrSTP)) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(timeEvaluated)) % Vectorsize == 0);

  static_assert(tensor::I<Cfg>::size() == tensor::Q<Cfg>::size(),
                "Sizes of tensors I and Q must match");

  kernel::derivativeTaylorExpansionEla<Cfg> krnl;
  krnl.I = timeEvaluated;
  const real* der = timeDerivativesOrSTP;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
    krnl.dQ(i) = der;
    der += tensor::dQ<Cfg>::size(i);
    krnl.power(i) = coeffs[i];
  }
  krnl.execute();
}

template <typename Cfg>
PerformanceEstimate Time<Cfg>::metrics() const {
  return PerformanceEstimate::fromKernel<kernel::derivativeTaylorExpansionEla<Cfg>>();
}

template <typename Cfg>
void Time<Cfg>::evaluateBatched(
    SEISSOL_GPU_PARAM const real* coeffs,
    SEISSOL_GPU_PARAM const real** timeDerivativesOrSTP,
    SEISSOL_GPU_PARAM real** timeIntegratedDofs,
    SEISSOL_GPU_PARAM std::size_t numElements,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  assert(timeDerivativesOrSTP != nullptr);
  assert(timeIntegratedDofs != nullptr);
  static_assert(tensor::I<Cfg>::size() == tensor::Q<Cfg>::size(),
                "Sizes of tensors I and Q must match");
  static_assert(kernel::gpu_derivativeTaylorExpansionEla<Cfg>::TmpMaxMemRequiredInBytes == 0);

  kernel::gpu_derivativeTaylorExpansionEla<Cfg> krnl;
  krnl.numElements = numElements;
  krnl.I = timeIntegratedDofs;
  std::size_t derivativeOffset = 0;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
    krnl.dQ(i) = timeDerivativesOrSTP;
    krnl.extraOffset_dQ(i) = derivativeOffset;
    derivativeOffset += tensor::dQ<Cfg>::size(i);
    krnl.power(i) = coeffs[i];
  }

  krnl.streamPtr = runtime.stream();
  krnl.execute();
#else
  logError() << "No GPU implementation provided";
#endif
}

template <typename Cfg>
void Spacetime<Cfg>::computeBatchedAder(
    SEISSOL_GPU_PARAM const real* coeffs,
    double /*timeStepWidth*/,
    LTS::Layer& /*layer*/,
    LocalTmp<Cfg>& /*tmp*/,
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& dataTable,
    bool /*updateDisplacement*/,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;
  /*
   * compute ADER scheme.
   */
  const ConditionalKey timeVolumeKernelKey(KernelNames::Time || KernelNames::Volume);
  if (dataTable.find(timeVolumeKernelKey) != dataTable.end()) {
    kernel::gpu_derivative<Cfg> krnl = deviceKrnlPrototype_;
    auto& entry = dataTable[timeVolumeKernelKey];

    const auto numElements = (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getSize();
    krnl.numElements = numElements;
    krnl.I = (entry.get<real*>(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr();
    krnl.Iane = (entry.get<real*>(inner_keys::Wp::Id::IdofsAne))->getDeviceDataPtr();

    std::size_t derivativesOffset = tensor::dQ<Cfg>::size(0);
    krnl.dQ(0) = (entry.get<real*>(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr();
    krnl.dQane(0) = (entry.get<real*>(inner_keys::Wp::Id::DofsAne))->getDeviceDataPtr();
    for (std::size_t i = 1; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
      krnl.dQ(i) = (entry.get<real*>(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr();
      krnl.extraOffset_dQ(i) = derivativesOffset;
      krnl.dQane(i) = (entry.get<real*>(inner_keys::Wp::Id::DerivativesAne))->getDeviceDataPtr();
      krnl.extraOffset_dQane(i) = i % 2 == 1 ? 0 : tensor::dQane<Cfg>::size(1);
      krnl.dQext(i) = (entry.get<real*>(inner_keys::Wp::Id::DerivativesExt))->getDeviceDataPtr();
      krnl.extraOffset_dQext(i) = i % 2 == 1 ? 0 : tensor::dQext<Cfg>::size(1);

      // TODO: compress
      derivativesOffset += tensor::dQ<Cfg>::size(i);
    }
    krnl.Q =
        const_cast<const real**>((entry.get<real*>(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr());

    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData<Cfg>, starMatrices);
    for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
      krnl.star(i) = const_cast<const real**>(
          (entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
      krnl.extraOffset_star(i) = SEISSOL_ARRAY_OFFSET(LocalIntegrationData<Cfg>, starMatrices, i);
    }

    krnl.W = const_cast<const real**>(
        entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    krnl.extraOffset_W = SEISSOL_OFFSET(LocalIntegrationData<Cfg>, specific.W);
    krnl.w = const_cast<const real**>(
        entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    krnl.extraOffset_w = SEISSOL_OFFSET(LocalIntegrationData<Cfg>, specific.w);
    krnl.E = const_cast<const real**>(
        entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    krnl.extraOffset_E = SEISSOL_OFFSET(LocalIntegrationData<Cfg>, specific.E);

    SEISSOL_OFFSET_ASSERT(LocalIntegrationData<Cfg>, specific.W);
    SEISSOL_OFFSET_ASSERT(LocalIntegrationData<Cfg>, specific.w);
    SEISSOL_OFFSET_ASSERT(LocalIntegrationData<Cfg>, specific.E);

    for (std::size_t der = 0; der < Cfg::ConvergenceOrder; ++der) {
      // update scalar for this derivative
      krnl.power(der) = coeffs[der];
    }

    krnl.streamPtr = runtime.stream();
    krnl.execute();
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

#define SEISSOL_INSTANTIATE(Cfg)                                                                   \
  template class Spacetime<Cfg>;                                                                   \
  template class Time<Cfg>;
SEISSOL_FOR_EACH_CONFIG_LINEARCKANELASTIC(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels::solver::linearckanelastic
