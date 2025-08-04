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

#include "Memory/GlobalData.h"

namespace seissol::kernels::solver::linearckanelastic {

template <typename Cfg>
void Local<Cfg>::setGlobalData(const GlobalData& global) {

#ifndef NDEBUG
  for (std::size_t stiffness = 0; stiffness < Cell::Dim; ++stiffness) {
    assert((reinterpret_cast<uintptr_t>(global.get<Cfg>().stiffnessMatrices(stiffness))) %
               Alignment ==
           0);
  }
  for (std::size_t flux = 0; flux < Cell::NumFaces; ++flux) {
    assert((reinterpret_cast<uintptr_t>(
               global.get<Cfg>().localChangeOfBasisMatricesTransposed(flux))) %
               Alignment ==
           0);
    assert((reinterpret_cast<uintptr_t>(global.get<Cfg>().changeOfBasisMatrices(flux))) %
               Alignment ==
           0);
  }
#endif

  m_volumeKernelPrototype.kDivM = global.get<Cfg>().stiffnessMatrices;
  m_localFluxKernelPrototype.rDivM = global.get<Cfg>().changeOfBasisMatrices;
  m_localFluxKernelPrototype.fMrT = global.get<Cfg>().localChangeOfBasisMatricesTransposed;
  m_localKernelPrototype.selectEla = init::selectEla<Cfg>::Values;
  m_localKernelPrototype.selectAne = init::selectAne<Cfg>::Values;

#ifdef ACL_DEVICE
  deviceVolumeKernelPrototype.kDivM = global.get<Cfg, Executor::Device>().stiffnessMatrices;
#ifdef USE_PREMULTIPLY_FLUX
  deviceLocalFluxKernelPrototype.plusFluxMatrices =
      global.get<Cfg, Executor::Device>().plusFluxMatrices;
#else
  deviceLocalFluxKernelPrototype.rDivM = global.get<Cfg, Executor::Device>().changeOfBasisMatrices;
  deviceLocalFluxKernelPrototype.fMrT =
      global.get<Cfg, Executor::Device>().localChangeOfBasisMatricesTransposed;
#endif
  deviceLocalKernelPrototype.selectEla = global.get<Cfg, Executor::Device>().selectEla;
  deviceLocalKernelPrototype.selectAne = global.get<Cfg, Executor::Device>().selectAne;
#endif
}

template <typename Cfg>
void Local<Cfg>::computeIntegral(real timeIntegratedDegreesOfFreedom[tensor::I<Cfg>::size()],
                                 LTS::Ref<Cfg>& data,
                                 LocalTmp<Cfg>& tmp,
                                 // TODO(Lukas) Nullable cause miniseissol. Maybe fix?
                                 const CellMaterialData* materialData,
                                 const CellBoundaryMapping<Cfg> (*cellBoundaryMapping)[4],
                                 double time,
                                 double timeStepWidth) {
  // assert alignments
#ifndef NDEBUG
  assert((reinterpret_cast<uintptr_t>(timeIntegratedDegreesOfFreedom)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(tmp.timeIntegratedAne)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(data.template get<LTS::Dofs>())) % Alignment == 0);
#endif

  alignas(Alignment) real Qext[tensor::Qext<Cfg>::size()];

  kernel::volumeExt<Cfg> volKrnl = m_volumeKernelPrototype;
  volKrnl.Qext = Qext;
  volKrnl.I = timeIntegratedDegreesOfFreedom;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
    volKrnl.star(i) = data.template get<LTS::LocalIntegration>().starMatrices[i];
  }

  kernel::localFluxExt<Cfg> lfKrnl = m_localFluxKernelPrototype;
  lfKrnl.Qext = Qext;
  lfKrnl.I = timeIntegratedDegreesOfFreedom;
  lfKrnl._prefetch.I = timeIntegratedDegreesOfFreedom + tensor::I<Cfg>::size();
  lfKrnl._prefetch.Q = data.template get<LTS::Dofs>() + tensor::Q<Cfg>::size();

  volKrnl.execute();

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // no element local contribution in the case of dynamic rupture boundary conditions
    if (data.template get<LTS::CellInformation>().faceTypes[face] != FaceType::DynamicRupture) {
      lfKrnl.AplusT = data.template get<LTS::LocalIntegration>().nApNm1[face];
      lfKrnl.execute(face);
    }
  }

  kernel::local<Cfg> lKrnl = m_localKernelPrototype;
  lKrnl.E = data.template get<LTS::LocalIntegration>().specific.E;
  lKrnl.Iane = tmp.timeIntegratedAne;
  lKrnl.Q = data.template get<LTS::Dofs>();
  lKrnl.Qane = data.template get<LTS::DofsAne>();
  lKrnl.Qext = Qext;
  lKrnl.W = data.template get<LTS::LocalIntegration>().specific.W;
  lKrnl.w = data.template get<LTS::LocalIntegration>().specific.w;

  lKrnl.execute();
}

template <typename Cfg>
void Local<Cfg>::flopsIntegral(const std::array<FaceType, Cell::NumFaces>& faceTypes,
                               std::uint64_t& nonZeroFlops,
                               std::uint64_t& hardwareFlops) {
  nonZeroFlops = seissol::kernel::volumeExt<Cfg>::NonZeroFlops;
  hardwareFlops = seissol::kernel::volumeExt<Cfg>::HardwareFlops;

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    if (faceTypes[face] != FaceType::DynamicRupture) {
      nonZeroFlops += seissol::kernel::localFluxExt<Cfg>::nonZeroFlops(face);
      hardwareFlops += seissol::kernel::localFluxExt<Cfg>::hardwareFlops(face);
    }
  }

  nonZeroFlops += seissol::kernel::local<Cfg>::NonZeroFlops;
  hardwareFlops += seissol::kernel::local<Cfg>::HardwareFlops;
}

template <typename Cfg>
std::uint64_t Local<Cfg>::bytesIntegral() {
  std::uint64_t reals = 0;

  // star matrices load
  reals += yateto::computeFamilySize<tensor::star<Cfg>>() + tensor::w<Cfg>::size() +
           tensor::W<Cfg>::size() + tensor::E<Cfg>::size();
  // flux solvers
  reals += 4 * tensor::AplusT<Cfg>::size();

  // DOFs write
  reals += tensor::Q<Cfg>::size() + tensor::Qane<Cfg>::size();

  return reals * sizeof(real);
}

template <typename Cfg>
void Local<Cfg>::computeBatchedIntegral(ConditionalPointersToRealsTable& dataTable,
                                        ConditionalMaterialTable& materialTable,
                                        ConditionalIndicesTable& indicesTable,
                                        double timeStepWidth,
                                        seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  // Volume integral
  ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  kernel::gpu_volumeExt<Cfg> volKrnl = deviceVolumeKernelPrototype;
  kernel::gpu_localFluxExt<Cfg> localFluxKrnl = deviceLocalFluxKernelPrototype;
  kernel::gpu_local<Cfg> localKrnl = deviceLocalKernelPrototype;

  if (dataTable.find(key) != dataTable.end()) {
    auto& entry = dataTable[key];

    volKrnl.numElements = (dataTable[key].get(inner_keys::Wp::Id::Dofs))->getSize();

    volKrnl.I =
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
    volKrnl.Qext = (entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();

    for (size_t i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
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

template <typename Cfg>
void Local<Cfg>::evaluateBatchedTimeDependentBc(
    ConditionalPointersToRealsTable& dataTable,
    ConditionalIndicesTable& indicesTable,
    LTS::Layer& layer,
    double time,
    double timeStepWidth,
    seissol::parallel::runtime::StreamRuntime& runtime) {}

#define _H_(cfg) template class Local<cfg>;
#include "ConfigInclude.h"

} // namespace seissol::kernels::solver::linearckanelastic
