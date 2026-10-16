// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "Local.h"

#include "Alignment.h"
#include "Common/Constants.h"
#include "Common/Marker.h"
#include "Config.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/Typedefs.h"
#include "Kernels/AnalyticalBoundary.h"
#include "Kernels/Interface.h"
#include "Memory/Descriptor/LTS.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

#include <array>
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
void Local<Cfg>::setGlobalData(const CompoundGlobalData<Cfg>& global) {
  volumeKernelPrototype_.bindGlobals(*global.onHost);
  localFluxKernelPrototype_.bindGlobals(*global.onHost);

  fsgFlux_.bindGlobals(*global.onHost);
  dirichletFlux_.bindGlobals(*global.onHost);
  nodalLfKrnlPrototype_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceVolumeKernelPrototype_.bindGlobals(*global.onDevice);
  deviceLocalFluxKernelPrototype_.bindGlobals(*global.onDevice);
  deviceFluxLocalAllKernelPrototype_.bindGlobals(*global.onDevice);
#endif
}

template <typename Cfg>
void Local<Cfg>::computeIntegral(real* timeIntegratedDoFs,
                                 LTS::Ref<Cfg>& data,
                                 LocalTmp<Cfg>& tmp,
                                 double time,
                                 double timeStepWidth) {
  // assert alignments
#ifndef NDEBUG
  assert((reinterpret_cast<uintptr_t>(timeIntegratedDoFs)) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(tmp.timeIntegratedAne)) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(data.template get<LTS::Dofs>())) % Vectorsize == 0);
#endif

  alignas(Alignment) real qext[tensor::Qext<Cfg>::size()];

  kernel::volumeExt<Cfg> volKrnl = volumeKernelPrototype_;
  volKrnl.Qext = qext;
  volKrnl.I = timeIntegratedDoFs;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
    volKrnl.star(i) = data.template get<LTS::LocalIntegration>().starMatrices[i];
  }

  kernel::localFluxExt<Cfg> lfKrnl = localFluxKernelPrototype_;
  lfKrnl.Qext = qext;
  lfKrnl.I = timeIntegratedDoFs;
  lfKrnl._prefetch.I = timeIntegratedDoFs + tensor::I<Cfg>::size();
  lfKrnl._prefetch.Q = data.template get<LTS::Dofs>() + tensor::Q<Cfg>::size();

  volKrnl.execute();

  const auto& cellBoundaryMapping = data.template get<LTS::BoundaryMapping>();
  const auto& materialData = data.template get<LTS::Material>();

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // no element local contribution in the case of dynamic rupture boundary conditions
    if (data.template get<LTS::CellInformation>().faceTypes[face] != FaceType::DynamicRupture) {
      lfKrnl.AplusT = data.template get<LTS::LocalIntegration>().nApNm1[face];
      lfKrnl.execute(face);
    }

    switch (data.template get<LTS::CellInformation>().faceTypes[face]) {
    case FaceType::FreeSurfaceGravity: {
      auto kernel = fsgFlux_;
      kernel.g2m = -2 * this->gravitationalAcceleration_;

      const real localRho = materialData.local->getDensity();
      kernel.rho = &localRho;
      kernel.averageNormalDisplacement = tmp.nodalAvgDisplacements[face].data();

      kernel.Qext = qext;
      kernel.AminusT = data.template get<LTS::NeighboringIntegration>().nAmNm1[face];

      kernel.execute(face);
      break;
    }
    case FaceType::Dirichlet: {
      auto kernel = dirichletFlux_;
      kernel.dirichletOffset = cellBoundaryMapping[face].dirichletOffset;
      kernel.dt = timeStepWidth;

      kernel.Qext = qext;
      kernel.AminusT = data.template get<LTS::NeighboringIntegration>().nAmNm1[face];

      kernel.execute(face);
      break;
    }
    case FaceType::Analytical: {
      assert(this->initConds_ != nullptr);
      const auto applyAnalyticalSolution =
          kernels::ApplyAnalyticalSolution<Cfg>(this->initConds_, data);

      alignas(Alignment) real dofsFaceBoundaryNodal[tensor::INodal<Cfg>::size()];
      analyticalBoundary_.evaluate(cellBoundaryMapping[face],
                                   applyAnalyticalSolution,
                                   dofsFaceBoundaryNodal,
                                   time,
                                   timeStepWidth);

      auto nodalLfKrnl = nodalLfKrnlPrototype_;
      nodalLfKrnl.Qext = qext;
      nodalLfKrnl.INodal = dofsFaceBoundaryNodal;
      nodalLfKrnl.AminusT = data.template get<LTS::NeighboringIntegration>().nAmNm1[face];
      nodalLfKrnl.execute(face);
      break;
    }
    default:
      // No boundary condition.
      break;
    }
  }

  kernel::local<Cfg> lKrnl = localKernelPrototype_;
  lKrnl.E = data.template get<LTS::LocalIntegration>().specific.E;
  lKrnl.Iane = tmp.timeIntegratedAne;
  lKrnl.Q = data.template get<LTS::Dofs>();
  lKrnl.Qane = data.template get<LTS::DofsAne>();
  lKrnl.Qext = qext;
  lKrnl.W = data.template get<LTS::LocalIntegration>().specific.W;
  lKrnl.w = data.template get<LTS::LocalIntegration>().specific.w;

  lKrnl.execute();
}

template <typename Cfg>
PerformanceEstimate
    Local<Cfg>::metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes) const {
  auto estimate = PerformanceEstimate::fromKernel<seissol::kernel::volumeExt<Cfg>>();

#if defined(ACL_DEVICE) && defined(SEISSOL_DEVICE_COMBINE_LOCAL_FLUX)
  constexpr bool CombineLocalFlux = true;
#else
  constexpr bool CombineLocalFlux = false;
#endif

  if constexpr (CombineLocalFlux) {
    estimate += PerformanceEstimate::fromKernel<seissol::kernel::fluxLocalAll<Cfg>>();
  } else {
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (faceTypes[face] != FaceType::DynamicRupture) {
        estimate += PerformanceEstimate::fromKernel<seissol::kernel::localFluxExt<Cfg>>(face);
      }
      switch (faceTypes[face]) {
      case FaceType::FreeSurfaceGravity:
        estimate += PerformanceEstimate::fromKernel<seissol::kernel::fsgFlux<Cfg>>(face);
        break;
      case FaceType::Dirichlet:
        estimate += PerformanceEstimate::fromKernel<seissol::kernel::dirichletFlux<Cfg>>(face);
        break;
      case FaceType::Analytical:
        estimate += PerformanceEstimate::fromKernel<seissol::kernel::localFluxNodal<Cfg>>(face);
        break;
      default:
        break;
      }
    }

    estimate += PerformanceEstimate::fromKernel<seissol::kernel::local<Cfg>>();
  }

  // legacy memory estimate
  std::uint64_t reals = 0;

  // star matrices load
  reals += yateto::computeFamilySize<tensor::star<Cfg>>() + tensor::w<Cfg>::size() +
           tensor::W<Cfg>::size() + tensor::E<Cfg>::size();
  // flux solvers
  reals += 4 * tensor::AplusT<Cfg>::size();

  // DOFs write
  reals += tensor::Q<Cfg>::size() + tensor::Qane<Cfg>::size();

  estimate.bytes = reals * sizeof(real);

  return estimate;
}

template <typename Cfg>
void Local<Cfg>::computeBatchedIntegral(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& dataTable,
    recording::ConditionalIndicesTable& /*indicesTable*/,
    double /*timeStepWidth*/,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;
  // Volume integral
  const ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  kernel::gpu_volumeExt<Cfg> volKrnl = deviceVolumeKernelPrototype_;

  if (dataTable.find(key) != dataTable.end()) {
    auto& entry = dataTable[key];

    volKrnl.numElements = (dataTable[key].get<real*>(inner_keys::Wp::Id::Dofs))->getSize();

    volKrnl.I =
        const_cast<const real**>((entry.get<real*>(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
    volKrnl.Qext = (entry.get<real*>(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();

    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData<Cfg>, starMatrices);
    for (size_t i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
      volKrnl.star(i) = const_cast<const real**>(
          (entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
      volKrnl.extraOffset_star(i) =
          SEISSOL_ARRAY_OFFSET(LocalIntegrationData<Cfg>, starMatrices, i);
    }
    volKrnl.streamPtr = runtime.stream();
    volKrnl.execute();

#ifdef SEISSOL_DEVICE_COMBINE_LOCAL_FLUX
    auto krnl = deviceFluxLocalAllKernelPrototype_;

    krnl.numElements = entry.get<real*>(inner_keys::Wp::Id::Dofs)->getSize();
    krnl.Q = (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    krnl.Qane = (entry.get<real*>(inner_keys::Wp::Id::DofsAne))->getDeviceDataPtr();
    krnl.Qext = (entry.get<real*>(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();
    krnl.Iane = const_cast<const real**>(
        (entry.get<real*>(inner_keys::Wp::Id::IdofsAne))->getDeviceDataPtr());
    krnl.W = const_cast<const real**>(
        entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    krnl.extraOffset_W = SEISSOL_OFFSET(LocalIntegrationData<Cfg>, specific.W);
    krnl.w = const_cast<const real**>(
        entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    krnl.extraOffset_w = SEISSOL_OFFSET(LocalIntegrationData<Cfg>, specific.w);
    krnl.E = const_cast<const real**>(
        entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    krnl.extraOffset_E = SEISSOL_OFFSET(LocalIntegrationData<Cfg>, specific.E);
    krnl.streamPtr = runtime.stream();

    SEISSOL_OFFSET_ASSERT(LocalIntegrationData<Cfg>, specific.W);
    SEISSOL_OFFSET_ASSERT(LocalIntegrationData<Cfg>, specific.w);
    SEISSOL_OFFSET_ASSERT(LocalIntegrationData<Cfg>, specific.E);

    krnl.I =
        const_cast<const real**>((entry.get<real*>(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());

    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData<Cfg>, nApNm1);
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      krnl.AplusTAll(face) = const_cast<const real**>(
          entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
      krnl.extraOffset_AplusTAll(face) =
          SEISSOL_ARRAY_OFFSET(LocalIntegrationData<Cfg>, nApNm1, face);
    }

    krnl.execute();
#endif
  }

// deprecated code; kept for comparison reasons
#ifndef SEISSOL_DEVICE_COMBINE_LOCAL_FLUX
  kernel::gpu_localFluxExt<Cfg> localFluxKrnl = deviceLocalFluxKernelPrototype_;
  kernel::gpu_local<Cfg> localKrnl = deviceLocalKernelPrototype_;

  // Local Flux Integral
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    key = ConditionalKey(*KernelNames::LocalFlux, !FaceKinds::DynamicRupture, face);

    if (dataTable.find(key) != dataTable.end()) {
      auto& entry = dataTable[key];
      localFluxKrnl.numElements = entry.get<real*>(inner_keys::Wp::Id::Dofs)->getSize();
      localFluxKrnl.Qext = (entry.get<real*>(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();
      localFluxKrnl.I = const_cast<const real**>(
          (entry.get<real*>(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
      localFluxKrnl.AplusT = const_cast<const real**>(
          entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());

      SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData<Cfg>, nApNm1);
      localFluxKrnl.extraOffset_AplusT =
          SEISSOL_ARRAY_OFFSET(LocalIntegrationData<Cfg>, nApNm1, face);
      localFluxKrnl.streamPtr = runtime.stream();
      localFluxKrnl.execute(face);
    }
  }

  key = ConditionalKey(KernelNames::Time || KernelNames::Volume);
  if (dataTable.find(key) != dataTable.end()) {
    auto& entry = dataTable[key];

    localKrnl.numElements = entry.get<real*>(inner_keys::Wp::Id::Dofs)->getSize();
    localKrnl.Q = (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    localKrnl.Qane = (entry.get<real*>(inner_keys::Wp::Id::DofsAne))->getDeviceDataPtr();
    localKrnl.Qext = const_cast<const real**>(
        (entry.get<real*>(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr());
    localKrnl.Iane = const_cast<const real**>(
        (entry.get<real*>(inner_keys::Wp::Id::IdofsAne))->getDeviceDataPtr());
    localKrnl.W = const_cast<const real**>(
        entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    localKrnl.extraOffset_W = SEISSOL_OFFSET(LocalIntegrationData<Cfg>, specific.W);
    localKrnl.w = const_cast<const real**>(
        entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    localKrnl.extraOffset_w = SEISSOL_OFFSET(LocalIntegrationData<Cfg>, specific.w);
    localKrnl.E = const_cast<const real**>(
        entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    localKrnl.extraOffset_E = SEISSOL_OFFSET(LocalIntegrationData<Cfg>, specific.E);
    localKrnl.streamPtr = runtime.stream();

    SEISSOL_OFFSET_ASSERT(LocalIntegrationData<Cfg>, specific.W);
    SEISSOL_OFFSET_ASSERT(LocalIntegrationData<Cfg>, specific.w);
    SEISSOL_OFFSET_ASSERT(LocalIntegrationData<Cfg>, specific.E);

    localKrnl.execute();
  }
#endif

#else
  logError() << "No GPU implementation provided";
#endif
}

template <typename Cfg>
void Local<Cfg>::evaluateBatchedTimeDependentBc(
    recording::ConditionalPointersToRealsTable& dataTable,
    recording::ConditionalIndicesTable& indicesTable,
    LTS::Layer& layer,
    double time,
    double timeStepWidth,
    seissol::parallel::runtime::StreamRuntime& runtime) {}

#define SEISSOL_INSTANTIATE(Cfg) template class Local<Cfg>;
SEISSOL_FOR_EACH_CONFIG_LINEARCKANELASTIC(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels::solver::linearckanelastic
