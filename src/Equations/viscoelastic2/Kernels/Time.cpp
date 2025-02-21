// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "Kernels/Time.h"

#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

#include <cassert>
#include <cstring>
#include <stdint.h>

#include <yateto.h>

#include "Kernels/DenseMatrixOps.h"
#include "generated_code/init.h"

namespace seissol::kernels {

void Time::setHostGlobalData(const GlobalData* global) {
  assert((reinterpret_cast<uintptr_t>(global->stiffnessMatricesTransposed(0))) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(global->stiffnessMatricesTransposed(1))) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(global->stiffnessMatricesTransposed(2))) % Alignment == 0);

  m_krnlPrototype.kDivMT = global->stiffnessMatricesTransposed;
  m_krnlPrototype.selectAne = init::selectAne::Values;
  m_krnlPrototype.selectEla = init::selectEla::Values;
}

void Time::setGlobalData(const CompoundGlobalData& global) {
  setHostGlobalData(global.onHost);

#ifdef ACL_DEVICE
  deviceKrnlPrototype.kDivMT = global.onDevice->stiffnessMatricesTransposed;
  // the selectAne/selectEla are inlined
  deviceKrnlPrototype.selectAne = global.onDevice->selectAne;
  deviceKrnlPrototype.selectEla = global.onDevice->selectEla;
#endif
}

void Time::computeAder(double timeStepWidth,
                       LocalData& data,
                       LocalTmp& tmp,
                       real timeIntegrated[tensor::I::size()],
                       real* timeDerivatives,
                       bool updateDisplacement) {
  /*
   * assert alignments.
   */
  assert((reinterpret_cast<uintptr_t>(data.dofs())) % Alignment == 0);
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

  krnl.dQ(0) = const_cast<real*>(data.dofs());
  if (timeDerivatives != nullptr) {
    streamstore(tensor::dQ::size(0), data.dofs(), timeDerivatives);
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

  krnl.dQane(0) = const_cast<real*>(data.dofsAne());
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQane(i) = temporaryBufferAne[i % 2];
    krnl.dQext(i) = temporaryBufferExt[i % 2];
  }

  krnl.I = timeIntegrated;
  krnl.Iane = tmp.timeIntegratedAne;

  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    krnl.star(i) = data.localIntegration().starMatrices[i];
  }
  krnl.w = data.localIntegration().specific.w;
  krnl.W = data.localIntegration().specific.W;
  krnl.E = data.localIntegration().specific.E;

  // powers in the taylor-series expansion
  krnl.power(0) = timeStepWidth;

  for (std::size_t der = 1; der < ConvergenceOrder; ++der) {
    // update scalar for this derivative
    krnl.power(der) = krnl.power(der - 1) * timeStepWidth / real(der + 1);
  }

  krnl.execute();

  // TODO(Lukas) Implement!
  // Compute integrated displacement over time step if needed.
}

void Time::flopsAder(unsigned int& nonZeroFlops, unsigned int& hardwareFlops) {
  nonZeroFlops = kernel::derivative::NonZeroFlops;
  hardwareFlops = kernel::derivative::HardwareFlops;
}

unsigned Time::bytesAder() {
  unsigned reals = 0;

  // DOFs load, tDOFs load, tDOFs write
  reals +=
      tensor::Q::size() + tensor::Qane::size() + 2 * tensor::I::size() + 2 * tensor::Iane::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star>() + tensor::w::size() + tensor::W::size() +
           tensor::E::size();

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

  // compute lengths of integration intervals
  real deltaTLower = integrationStart - expansionPoint;
  real deltaTUpper = integrationEnd - expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  real firstTerm = static_cast<real>(1.0);
  real secondTerm = static_cast<real>(1.0);
  real factorial = static_cast<real>(1.0);

  kernel::derivativeTaylorExpansionEla intKrnl;
  intKrnl.I = timeIntegrated;
  const real* der = timeDerivatives;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = der;
    der += tensor::dQ::size(i);
  }

  // iterate over time derivatives
  for (std::size_t der = 0; der < ConvergenceOrder; ++der) {
    firstTerm *= deltaTUpper;
    secondTerm *= deltaTLower;
    factorial *= static_cast<real>(der + 1);

    intKrnl.power(der) = firstTerm - secondTerm;
    intKrnl.power(der) /= factorial;
  }

  intKrnl.execute();
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

  real deltaT = time - expansionPoint;

  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");

  kernel::derivativeTaylorExpansionEla intKrnl;
  intKrnl.I = timeEvaluated;
  const real* der = timeDerivatives;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = der;
    der += tensor::dQ::size(i);
  }
  intKrnl.power(0) = 1.0;

  // iterate over time derivatives
  for (std::size_t derivative = 1; derivative < ConvergenceOrder; ++derivative) {
    intKrnl.power(derivative) = intKrnl.power(derivative - 1) * deltaT / real(derivative);
  }

  intKrnl.execute();
}

void Time::flopsTaylorExpansion(long long& nonZeroFlops, long long& hardwareFlops) {
  // interate over derivatives
  nonZeroFlops = kernel::derivativeTaylorExpansionEla::NonZeroFlops;
  hardwareFlops = kernel::derivativeTaylorExpansionEla::HardwareFlops;
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
  real deltaTLower = integrationStart - expansionPoint;
  real deltaTUpper = integrationEnd - expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  real firstTerm = static_cast<real>(1.0);
  real secondTerm = static_cast<real>(1.0);
  real factorial = static_cast<real>(1.0);

  kernel::gpu_derivativeTaylorExpansionEla intKrnl;
  intKrnl.numElements = numElements;
  real* tmpMem = reinterpret_cast<real*>(device.api->getStackMemory(
      seissol::kernel::gpu_derivativeTaylorExpansionEla::TmpMaxMemRequiredInBytes * numElements));

  intKrnl.I = timeIntegratedDofs;

  unsigned derivativesOffset = 0;
  for (size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = timeDerivatives;
    intKrnl.extraOffset_dQ(i) = derivativesOffset;
    derivativesOffset += tensor::dQ::size(i);
  }

  // iterate over time derivatives
  for (int der = 0; der < ConvergenceOrder; ++der) {
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
  assert(false && "no implementation provided");
#endif
}

void seissol::kernels::Time::computeBatchedTaylorExpansion(
    real time,
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
  static_assert(kernel::gpu_derivativeTaylorExpansionEla::TmpMaxMemRequiredInBytes == 0);

  kernel::gpu_derivativeTaylorExpansionEla intKrnl;
  intKrnl.numElements = numElements;
  intKrnl.I = timeEvaluated;
  std::size_t derivativeOffset = 0;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = const_cast<const real**>(timeDerivatives);
    intKrnl.extraOffset_dQ(i) = derivativeOffset;
    derivativeOffset += tensor::dQ::size(i);
  }

  // iterate over time derivatives
  const real deltaT = time - expansionPoint;
  intKrnl.power(0) = 1.0;
  for (int derivative = 1; derivative < ConvergenceOrder; ++derivative) {
    intKrnl.power(derivative) =
        intKrnl.power(derivative - 1) * deltaT / static_cast<real>(derivative);
  }

  intKrnl.streamPtr = runtime.stream();
  intKrnl.execute();
#else
  assert(false && "no implementation provided");
#endif
}

void Time::computeBatchedAder(double timeStepWidth,
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

    unsigned starOffset = 0;
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      krnl.star(i) =
          const_cast<const real**>((entry.get(inner_keys::Wp::Id::Star))->getDeviceDataPtr());
      krnl.extraOffset_star(i) = starOffset;
      starOffset += tensor::star::size(i);
    }

    krnl.w = const_cast<const real**>((entry.get(inner_keys::Wp::Id::Omega))->getDeviceDataPtr());
    krnl.W = const_cast<const real**>((entry.get(inner_keys::Wp::Id::W))->getDeviceDataPtr());
    krnl.E = const_cast<const real**>((entry.get(inner_keys::Wp::Id::E))->getDeviceDataPtr());

    // powers in the taylor-series expansion
    krnl.power(0) = timeStepWidth;

    for (unsigned der = 1; der < ConvergenceOrder; ++der) {
      // update scalar for this derivative
      krnl.power(der) = krnl.power(der - 1) * timeStepWidth / real(der + 1);
    }

    device.algorithms.streamBatchedData(
        (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr(),
        (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr(),
        tensor::Q::Size,
        krnl.numElements,
        runtime.stream());

    krnl.streamPtr = runtime.stream();
    krnl.execute();
  }
#else
  assert(false && "no implementation provided");
#endif
}

} // namespace seissol::kernels
