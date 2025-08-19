// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "TimeBase.h"

#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

#include <cassert>
#include <cstring>
#include <stdint.h>

#include <yateto.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#endif

#include "Kernels/MemoryOps.h"
#include "generated_code/init.h"

namespace seissol::kernels::solver::linearckanelastic {

void Time::setGlobalData(const CompoundGlobalData& global) {}

void Spacetime::setGlobalData(const CompoundGlobalData& global) {
  assert((reinterpret_cast<uintptr_t>(global.onHost->stiffnessMatricesTransposed(0))) % Alignment ==
         0);
  assert((reinterpret_cast<uintptr_t>(global.onHost->stiffnessMatricesTransposed(1))) % Alignment ==
         0);
  assert((reinterpret_cast<uintptr_t>(global.onHost->stiffnessMatricesTransposed(2))) % Alignment ==
         0);

  m_krnlPrototype.kDivMT = global.onHost->stiffnessMatricesTransposed;
  m_krnlPrototype.selectAne = init::selectAne::Values;
  m_krnlPrototype.selectEla = init::selectEla::Values;

#ifdef ACL_DEVICE
  deviceKrnlPrototype.kDivMT = global.onDevice->stiffnessMatricesTransposed;
  // the selectAne/selectEla are inlined
  deviceKrnlPrototype.selectAne = global.onDevice->selectAne;
  deviceKrnlPrototype.selectEla = global.onDevice->selectEla;
#endif
}

void Spacetime::computeAder(const real* coeffs,
                            double timeStepWidth,
                            LTS::Ref& data,
                            LocalTmp& tmp,
                            real timeIntegrated[tensor::I::size()],
                            real* timeDerivatives,
                            bool updateDisplacement) {
  /*
   * assert alignments.
   */
  assert((reinterpret_cast<uintptr_t>(data.get<LTS::Dofs>())) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(timeIntegrated)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(timeDerivatives)) % Alignment == 0 ||
         timeDerivatives == NULL);

  /*
   * compute ADER scheme.
   */
  // temporary result
  // TODO(David): move these temporary buffers into the Yateto kernel, maybe
  alignas(PagesizeStack) real temporaryBuffer[2][tensor::dQ::size(0)];
  alignas(PagesizeStack) real temporaryBufferExt[2][tensor::dQext::size(1)];
  alignas(PagesizeStack) real temporaryBufferAne[2][tensor::dQane::size(0)];

  kernel::derivative krnl = m_krnlPrototype;

  krnl.dQ(0) = const_cast<real*>(data.get<LTS::Dofs>());
  if (timeDerivatives != nullptr) {
    streamstore(tensor::dQ::size(0), data.get<LTS::Dofs>(), timeDerivatives);
    real* derOut = timeDerivatives;
    for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
      derOut += tensor::dQ::size(i - 1);
      krnl.dQ(i) = derOut;
    }
  } else {
    for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
      krnl.dQ(i) = temporaryBuffer[i % 2];
    }
  }

  krnl.dQane(0) = const_cast<real*>(data.get<LTS::DofsAne>());
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQane(i) = temporaryBufferAne[i % 2];
    krnl.dQext(i) = temporaryBufferExt[i % 2];
  }

  krnl.I = timeIntegrated;
  krnl.Iane = tmp.timeIntegratedAne;

  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    krnl.star(i) = data.get<LTS::LocalIntegration>().starMatrices[i];
  }
  krnl.w = data.get<LTS::LocalIntegration>().specific.w;
  krnl.W = data.get<LTS::LocalIntegration>().specific.W;
  krnl.E = data.get<LTS::LocalIntegration>().specific.E;

  // powers in the taylor-series expansion
  for (std::size_t der = 0; der < ConvergenceOrder; ++der) {
    krnl.power(der) = coeffs[der];
  }

  krnl.execute();

  // TODO(Lukas) Implement!
  // Compute integrated displacement over time step if needed.
}

void Spacetime::flopsAder(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) {
  nonZeroFlops = kernel::derivative::NonZeroFlops;
  hardwareFlops = kernel::derivative::HardwareFlops;
}

std::uint64_t Spacetime::bytesAder() {
  std::uint64_t reals = 0;

  // DOFs load, tDOFs load, tDOFs write
  reals +=
      tensor::Q::size() + tensor::Qane::size() + 2 * tensor::I::size() + 2 * tensor::Iane::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star>() + tensor::w::size() + tensor::W::size() +
           tensor::E::size();

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

  kernel::derivativeTaylorExpansionEla krnl;
  krnl.I = timeEvaluated;
  const real* der = timeDerivatives;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQ(i) = der;
    der += tensor::dQ::size(i);
    krnl.power(i) = coeffs[i];
  }
  krnl.execute();
}

void Time::flopsEvaluate(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) {
  // interate over derivatives
  nonZeroFlops = kernel::derivativeTaylorExpansionEla::NonZeroFlops;
  hardwareFlops = kernel::derivativeTaylorExpansionEla::HardwareFlops;
}

void Time::evaluateBatched(const real* coeffs,
                           const real** timeDerivatives,
                           real** timeIntegratedDofs,
                           std::size_t numElements,
                           seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  assert(timeDerivatives != nullptr);
  assert(timeIntegratedDofs != nullptr);
  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");
  static_assert(kernel::gpu_derivativeTaylorExpansionEla::TmpMaxMemRequiredInBytes == 0);

  kernel::gpu_derivativeTaylorExpansionEla krnl;
  krnl.numElements = numElements;
  krnl.I = timeIntegratedDofs;
  std::size_t derivativeOffset = 0;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQ(i) = const_cast<const real**>(timeDerivatives);
    krnl.extraOffset_dQ(i) = derivativeOffset;
    derivativeOffset += tensor::dQ::size(i);
    krnl.power(i) = coeffs[i];
  }

  krnl.streamPtr = runtime.stream();
  krnl.execute();
#else
  logError() << "No GPU implementation provided";
#endif
}

void Spacetime::computeBatchedAder(const real* coeffs,
                                   double timeStepWidth,
                                   LocalTmp& tmp,
                                   ConditionalPointersToRealsTable& dataTable,
                                   ConditionalMaterialTable& materialTable,
                                   bool updateDisplacement,
                                   seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  /*
   * compute ADER scheme.
   */
  ConditionalKey timeVolumeKernelKey(KernelNames::Time || KernelNames::Volume);
  if (dataTable.find(timeVolumeKernelKey) != dataTable.end()) {
    kernel::gpu_derivative krnl = deviceKrnlPrototype;
    auto& entry = dataTable[timeVolumeKernelKey];

    const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
    krnl.numElements = numElements;
    krnl.I = (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr();
    krnl.Iane = (entry.get(inner_keys::Wp::Id::IdofsAne))->getDeviceDataPtr();

    unsigned derivativesOffset = tensor::dQ::size(0);
    krnl.dQ(0) = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    krnl.dQane(0) = (entry.get(inner_keys::Wp::Id::DofsAne))->getDeviceDataPtr();
    for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
      krnl.dQ(i) = (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr();
      krnl.extraOffset_dQ(i) = derivativesOffset;
      krnl.dQane(i) = (entry.get(inner_keys::Wp::Id::DerivativesAne))->getDeviceDataPtr();
      krnl.extraOffset_dQane(i) = i % 2 == 1 ? 0 : tensor::dQane::size(1);
      krnl.dQext(i) = (entry.get(inner_keys::Wp::Id::DerivativesExt))->getDeviceDataPtr();
      krnl.extraOffset_dQext(i) = i % 2 == 1 ? 0 : tensor::dQext::size(1);

      // TODO: compress
      derivativesOffset += tensor::dQ::size(i);
    }

    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      krnl.star(i) = const_cast<const real**>(
          (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
      krnl.extraOffset_star(i) = SEISSOL_ARRAY_OFFSET(LocalIntegrationData, starMatrices, i);
    }

    krnl.W = const_cast<const real**>(
        entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    krnl.extraOffset_W = SEISSOL_OFFSET(LocalIntegrationData, specific.W);
    krnl.w = const_cast<const real**>(
        entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    krnl.extraOffset_w = SEISSOL_OFFSET(LocalIntegrationData, specific.w);
    krnl.E = const_cast<const real**>(
        entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    krnl.extraOffset_E = SEISSOL_OFFSET(LocalIntegrationData, specific.E);

    for (std::size_t der = 0; der < ConvergenceOrder; ++der) {
      // update scalar for this derivative
      krnl.power(der) = coeffs[der];
    }

    device.algorithms.streamBatchedData(
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr()),
        (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr(),
        tensor::Q::Size,
        krnl.numElements,
        runtime.stream());

    krnl.streamPtr = runtime.stream();
    krnl.execute();
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

} // namespace seissol::kernels::solver::linearckanelastic
