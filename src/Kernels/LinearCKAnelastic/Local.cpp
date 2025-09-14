// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "LocalBase.h"

#include <cassert>
#include <cstring>
#include <stdint.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#endif

#include <yateto.h>

namespace seissol::kernels::solver::linearckanelastic {

void Local::setGlobalData(const CompoundGlobalData& global) {

#ifndef NDEBUG
  for (std::size_t stiffness = 0; stiffness < Cell::Dim; ++stiffness) {
    assert((reinterpret_cast<uintptr_t>(global.onHost->stiffnessMatrices(stiffness))) % Alignment ==
           0);
  }
  for (std::size_t flux = 0; flux < Cell::NumFaces; ++flux) {
    assert(
        (reinterpret_cast<uintptr_t>(global.onHost->localChangeOfBasisMatricesTransposed(flux))) %
            Alignment ==
        0);
    assert((reinterpret_cast<uintptr_t>(global.onHost->changeOfBasisMatrices(flux))) % Alignment ==
           0);
  }
#endif

  m_volumeKernelPrototype.kDivM = global.onHost->stiffnessMatrices;
  m_localFluxKernelPrototype.rDivM = global.onHost->changeOfBasisMatrices;
  m_localFluxKernelPrototype.fMrT = global.onHost->localChangeOfBasisMatricesTransposed;

#ifdef ACL_DEVICE
  deviceVolumeKernelPrototype.kDivM = global.onDevice->stiffnessMatrices;
#ifdef USE_PREMULTIPLY_FLUX
  deviceLocalFluxKernelPrototype.plusFluxMatrices = global.onDevice->plusFluxMatrices;
#else
  deviceLocalFluxKernelPrototype.rDivM = global.onDevice->changeOfBasisMatrices;
  deviceLocalFluxKernelPrototype.fMrT = global.onDevice->localChangeOfBasisMatricesTransposed;
#endif
#endif
}

void Local::computeIntegral(real timeIntegratedDegreesOfFreedom[tensor::I::size()],
                            LocalData& data,
                            LocalTmp& tmp,
                            // TODO(Lukas) Nullable cause miniseissol. Maybe fix?
                            const CellMaterialData* materialData,
                            const CellBoundaryMapping (*cellBoundaryMapping)[4],
                            double time,
                            double timeStepWidth) {
  // assert alignments
#ifndef NDEBUG
  assert((reinterpret_cast<uintptr_t>(timeIntegratedDegreesOfFreedom)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(tmp.timeIntegratedAne)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(data.dofs())) % Alignment == 0);
#endif

  alignas(Alignment) real Qext[tensor::Qext::size()];

  kernel::volumeExt volKrnl = m_volumeKernelPrototype;
  volKrnl.Qext = Qext;
  volKrnl.I = timeIntegratedDegreesOfFreedom;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    volKrnl.star(i) = data.localIntegration().starMatrices[i];
  }

  kernel::localFluxExt lfKrnl = m_localFluxKernelPrototype;
  lfKrnl.Qext = Qext;
  lfKrnl.I = timeIntegratedDegreesOfFreedom;
  lfKrnl._prefetch.I = timeIntegratedDegreesOfFreedom + tensor::I::size();
  lfKrnl._prefetch.Q = data.dofs() + tensor::Q::size();

  volKrnl.execute();

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // no element local contribution in the case of dynamic rupture boundary conditions
    if (data.cellInformation().faceTypes[face] != FaceType::DynamicRupture) {
      lfKrnl.AplusT = data.localIntegration().nApNm1[face];
      lfKrnl.execute(face);
    }
  }

  kernel::local lKrnl = m_localKernelPrototype;
  lKrnl.E = data.localIntegration().specific.E;
  lKrnl.Iane = tmp.timeIntegratedAne;
  lKrnl.Q = data.dofs();
  lKrnl.Qane = data.dofsAne();
  lKrnl.Qext = Qext;
  lKrnl.W = data.localIntegration().specific.W;
  lKrnl.w = data.localIntegration().specific.w;

  lKrnl.execute();
}

void Local::flopsIntegral(const FaceType faceTypes[4],
                          std::uint64_t& nonZeroFlops,
                          std::uint64_t& hardwareFlops) {
  nonZeroFlops = seissol::kernel::volumeExt::NonZeroFlops;
  hardwareFlops = seissol::kernel::volumeExt::HardwareFlops;

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    if (faceTypes[face] != FaceType::DynamicRupture) {
      nonZeroFlops += seissol::kernel::localFluxExt::nonZeroFlops(face);
      hardwareFlops += seissol::kernel::localFluxExt::hardwareFlops(face);
    }
  }

  nonZeroFlops += seissol::kernel::local::NonZeroFlops;
  hardwareFlops += seissol::kernel::local::HardwareFlops;
}

std::uint64_t Local::bytesIntegral() {
  std::uint64_t reals = 0;

  // star matrices load
  reals += yateto::computeFamilySize<tensor::star>() + tensor::w::size() + tensor::W::size() +
           tensor::E::size();
  // flux solvers
  reals += 4 * tensor::AplusT::size();

  // DOFs write
  reals += tensor::Q::size() + tensor::Qane::size();

  return reals * sizeof(real);
}

void Local::computeBatchedIntegral(ConditionalPointersToRealsTable& dataTable,
                                   ConditionalMaterialTable& materialTable,
                                   ConditionalIndicesTable& indicesTable,
                                   kernels::LocalData::Loader& loader,
                                   LocalTmp& tmp,
                                   double timeStepWidth,
                                   seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  // Volume integral
  ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  kernel::gpu_volumeExt volKrnl = deviceVolumeKernelPrototype;
  kernel::gpu_localFluxExt localFluxKrnl = deviceLocalFluxKernelPrototype;
  kernel::gpu_local localKrnl = deviceLocalKernelPrototype;

  if (dataTable.find(key) != dataTable.end()) {
    auto& entry = dataTable[key];

    volKrnl.numElements = (dataTable[key].get(inner_keys::Wp::Id::Dofs))->getSize();

    volKrnl.I =
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
    volKrnl.Qext = (entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();

    for (size_t i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      volKrnl.star(i) = const_cast<const real**>(
          (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
      volKrnl.extraOffset_star(i) = SEISSOL_ARRAY_OFFSET(LocalIntegrationData, starMatrices, i);
    }
    volKrnl.streamPtr = runtime.stream();
    volKrnl.execute();
  }

  // Local Flux Integral
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    key = ConditionalKey(*KernelNames::LocalFlux, !FaceKinds::DynamicRupture, face);

    if (dataTable.find(key) != dataTable.end()) {
      auto& entry = dataTable[key];
      localFluxKrnl.numElements = entry.get(inner_keys::Wp::Id::Dofs)->getSize();
      localFluxKrnl.Qext = (entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();
      localFluxKrnl.I =
          const_cast<const real**>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
      localFluxKrnl.AplusT = const_cast<const real**>(
          entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
      localFluxKrnl.extraOffset_AplusT = SEISSOL_ARRAY_OFFSET(LocalIntegrationData, nApNm1, face);
      localFluxKrnl.streamPtr = runtime.stream();
      localFluxKrnl.execute(face);
    }
  }

  key = ConditionalKey(KernelNames::Time || KernelNames::Volume);
  if (dataTable.find(key) != dataTable.end()) {
    auto& entry = dataTable[key];

    localKrnl.numElements = entry.get(inner_keys::Wp::Id::Dofs)->getSize();
    localKrnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    localKrnl.Qane = (entry.get(inner_keys::Wp::Id::DofsAne))->getDeviceDataPtr();
    localKrnl.Qext =
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr());
    localKrnl.Iane =
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::IdofsAne))->getDeviceDataPtr());
    localKrnl.W = const_cast<const real**>(
        entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    localKrnl.extraOffset_W = SEISSOL_OFFSET(LocalIntegrationData, specific.W);
    localKrnl.w = const_cast<const real**>(
        entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    localKrnl.extraOffset_w = SEISSOL_OFFSET(LocalIntegrationData, specific.w);
    localKrnl.E = const_cast<const real**>(
        entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    localKrnl.extraOffset_E = SEISSOL_OFFSET(LocalIntegrationData, specific.E);
    localKrnl.streamPtr = runtime.stream();
    localKrnl.execute();
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

void Local::evaluateBatchedTimeDependentBc(ConditionalPointersToRealsTable& dataTable,
                                           ConditionalIndicesTable& indicesTable,
                                           kernels::LocalData::Loader& loader,
                                           seissol::initializer::Layer& layer,
                                           seissol::initializer::LTS& lts,
                                           double time,
                                           double timeStepWidth,
                                           seissol::parallel::runtime::StreamRuntime& runtime) {}

} // namespace seissol::kernels::solver::linearckanelastic
