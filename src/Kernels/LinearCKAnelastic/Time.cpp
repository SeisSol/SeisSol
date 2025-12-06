// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "Time.h"

#include "Common/Marker.h"
#include "GeneratedCode/init.h"
#include "Kernels/MemoryOps.h"

#include <Alignment.h>
#include <GeneratedCode/metagen/kernel.h>
#include <GeneratedCode/metagen/tensor.h>
#include <Kernels/Interface.h>
#include <Memory/Descriptor/LTS.h>
#include <Parallel/Runtime/Stream.h>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <stdint.h>
#include <utils/logger.h>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#endif

#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

#include "Memory/GlobalData.h"

namespace seissol::kernels::solver::linearckanelastic {

template <typename Cfg>
void Time<Cfg>::setGlobalData(const GlobalData& global) {}

template <typename Cfg>
void Spacetime<Cfg>::setGlobalData(const GlobalData& global) {
  assert((reinterpret_cast<uintptr_t>(global.get<Cfg>().stiffnessMatricesTransposed(0))) %
             Alignment ==
         0);
  assert((reinterpret_cast<uintptr_t>(global.get<Cfg>().stiffnessMatricesTransposed(1))) %
             Alignment ==
         0);
  assert((reinterpret_cast<uintptr_t>(global.get<Cfg>().stiffnessMatricesTransposed(2))) %
             Alignment ==
         0);

  m_krnlPrototype.kDivMT = global.get<Cfg>().stiffnessMatricesTransposed;
  m_krnlPrototype.selectAne = init::selectAne<Cfg>::Values;
  m_krnlPrototype.selectEla = init::selectEla<Cfg>::Values;

#ifdef ACL_DEVICE
  deviceKrnlPrototype.kDivMT = global.get<Cfg, Executor::Device>().stiffnessMatricesTransposed;
  // the selectAne/selectEla are inlined
  deviceKrnlPrototype.selectAne = global.get<Cfg, Executor::Device>().selectAne;
  deviceKrnlPrototype.selectEla = global.get<Cfg, Executor::Device>().selectEla;
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
  /*
   * assert alignments.
   */
  assert((reinterpret_cast<uintptr_t>(data.template get<LTS::Dofs>())) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(timeIntegrated)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(timeDerivatives)) % Alignment == 0 ||
         timeDerivatives == NULL);

  /*
   * compute ADER scheme.
   */
  // temporary result
  // TODO(David): move these temporary buffers into the Yateto kernel, maybe
  alignas(PagesizeStack) real temporaryBuffer[2][tensor::dQ<Cfg>::size(0)];
  alignas(PagesizeStack) real temporaryBufferExt[2][tensor::dQext<Cfg>::size(1)];
  alignas(PagesizeStack) real temporaryBufferAne[2][tensor::dQane<Cfg>::size(0)];

  kernel::derivative<Cfg> krnl = m_krnlPrototype;

  krnl.dQ(0) = const_cast<real*>(data.template get<LTS::Dofs>());
  if (timeDerivatives != nullptr) {
    streamstore(tensor::dQ<Cfg>::size(0), data.template get<LTS::Dofs>(), timeDerivatives);
    real* derOut = timeDerivatives;
    for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
      derOut += tensor::dQ<Cfg>::size(i - 1);
      krnl.dQ(i) = derOut;
    }
  } else {
    for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
      krnl.dQ(i) = temporaryBuffer[i % 2];
    }
  }

  krnl.dQane(0) = const_cast<real*>(data.template get<LTS::DofsAne>());
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
    krnl.dQane(i) = temporaryBufferAne[i % 2];
    krnl.dQext(i) = temporaryBufferExt[i % 2];
  }

  krnl.I = timeIntegrated;
  krnl.Iane = tmp.timeIntegratedAne;

  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
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

  // TODO(Lukas) Implement!
  // Compute integrated displacement over time step if needed.
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
  reals += tensor::Q<Cfg>::size() + tensor::Qane<Cfg>::size() + 2 * tensor::I<Cfg>::size() +
           2 * tensor::Iane<Cfg>::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star<Cfg>>() + tensor::w<Cfg>::size() +
           tensor::W<Cfg>::size() + tensor::E<Cfg>::size();

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

  kernel::derivativeTaylorExpansionEla<Cfg> krnl;
  krnl.I = timeEvaluated;
  const real* der = timeDerivatives;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
    krnl.dQ(i) = der;
    der += tensor::dQ<Cfg>::size(i);
    krnl.power(i) = coeffs[i];
  }
  krnl.execute();
}

template <typename Cfg>
void Time<Cfg>::flopsEvaluate(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) {
  // interate over derivatives
  nonZeroFlops = kernel::derivativeTaylorExpansionEla<Cfg>::NonZeroFlops;
  hardwareFlops = kernel::derivativeTaylorExpansionEla<Cfg>::HardwareFlops;
}

template <typename Cfg>
void Time<Cfg>::evaluateBatched(
    SEISSOL_GPU_PARAM const real* coeffs,
    SEISSOL_GPU_PARAM const real** timeDerivatives,
    SEISSOL_GPU_PARAM real** timeIntegratedDofs,
    SEISSOL_GPU_PARAM std::size_t numElements,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  assert(timeDerivatives != nullptr);
  assert(timeIntegratedDofs != nullptr);
  static_assert(tensor::I<Cfg>::size() == tensor::Q<Cfg>::size(),
                "Sizes of tensors I and Q must match");
  static_assert(kernel::gpu_derivativeTaylorExpansionEla<Cfg>::TmpMaxMemRequiredInBytes == 0);

  kernel::gpu_derivativeTaylorExpansionEla<Cfg> krnl;
  krnl.numElements = numElements;
  krnl.I = timeIntegratedDofs;
  std::size_t derivativeOffset = 0;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
    krnl.dQ(i) = const_cast<const real**>(timeDerivatives);
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
    SEISSOL_GPU_PARAM double timeStepWidth,
    SEISSOL_GPU_PARAM LocalTmp<Cfg>& tmp,
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& dataTable,
    SEISSOL_GPU_PARAM recording::ConditionalMaterialTable& materialTable,
    SEISSOL_GPU_PARAM bool updateDisplacement,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;
  /*
   * compute ADER scheme.
   */
  ConditionalKey timeVolumeKernelKey(KernelNames::Time || KernelNames::Volume);
  if (dataTable.find(timeVolumeKernelKey) != dataTable.end()) {
    kernel::gpu_derivative<Cfg> krnl = deviceKrnlPrototype;
    auto& entry = dataTable[timeVolumeKernelKey];

    const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
    krnl.numElements = numElements;
    krnl.I = (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtrAs<real*>();
    krnl.Iane = (entry.get(inner_keys::Wp::Id::IdofsAne))->getDeviceDataPtrAs<real*>();

    unsigned derivativesOffset = tensor::dQ<Cfg>::size(0);
    krnl.dQ(0) = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtrAs<real*>();
    krnl.dQane(0) = (entry.get(inner_keys::Wp::Id::DofsAne))->getDeviceDataPtrAs<real*>();
    for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
      krnl.dQ(i) = (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtrAs<real*>();
      krnl.extraOffset_dQ(i) = derivativesOffset;
      krnl.dQane(i) = (entry.get(inner_keys::Wp::Id::DerivativesAne))->getDeviceDataPtrAs<real*>();
      krnl.extraOffset_dQane(i) = i % 2 == 1 ? 0 : tensor::dQane<Cfg>::size(1);
      krnl.dQext(i) = (entry.get(inner_keys::Wp::Id::DerivativesExt))->getDeviceDataPtrAs<real*>();
      krnl.extraOffset_dQext(i) = i % 2 == 1 ? 0 : tensor::dQext<Cfg>::size(1);

      // TODO: compress
      derivativesOffset += tensor::dQ<Cfg>::size(i);
    }

    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
      krnl.star(i) = const_cast<const real**>(
          (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtrAs<real*>());
      krnl.extraOffset_star(i) = SEISSOL_ARRAY_OFFSET(LocalIntegrationData<Cfg>, starMatrices, i);
    }

    krnl.W = const_cast<const real**>(
        entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtrAs<real*>());
    krnl.extraOffset_W = SEISSOL_OFFSET(LocalIntegrationData<Cfg>, specific.W);
    krnl.w = const_cast<const real**>(
        entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtrAs<real*>());
    krnl.extraOffset_w = SEISSOL_OFFSET(LocalIntegrationData<Cfg>, specific.w);
    krnl.E = const_cast<const real**>(
        entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtrAs<real*>());
    krnl.extraOffset_E = SEISSOL_OFFSET(LocalIntegrationData<Cfg>, specific.E);

    for (std::size_t der = 0; der < Cfg::ConvergenceOrder; ++der) {
      // update scalar for this derivative
      krnl.power(der) = coeffs[der];
    }

    device.algorithms.streamBatchedData(
        const_cast<const real**>(
            (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtrAs<real*>()),
        (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtrAs<real*>(),
        tensor::Q<Cfg>::Size,
        krnl.numElements,
        runtime.stream());

    krnl.streamPtr = runtime.stream();
    krnl.execute();
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

#define SEISSOL_CONFIGITER(cfg) template class Spacetime<cfg>;
#include "ConfigIncludeLinearCKAne.h"

#define SEISSOL_CONFIGITER(cfg) template class Time<cfg>;
#include "ConfigIncludeLinearCKAne.h"

} // namespace seissol::kernels::solver::linearckanelastic
