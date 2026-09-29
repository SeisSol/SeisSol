// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "Local.h"

#include "Common/Marker.h"
#include "Initializer/Typedefs.h"
#include "Kernels/AnalyticalBoundary.h"
#include "Kernels/Common.h"
#include "Kernels/StarOperands.h"
#include "Monitoring/Metric.h"

#include <cassert>
#include <cstddef>
#include <cstring>
#include <stdint.h>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#endif

GENERATE_HAS_MEMBER(E)
GENERATE_HAS_MEMBER(extraOffset_E)

namespace seissol::kernels::solver::linearckanelastic {

void Local::setGlobalData(const CompoundGlobalData& global) {
  volumeKernelPrototype_.bindGlobals(*global.onHost);
  localFluxKernelPrototype_.bindGlobals(*global.onHost);
  // the relaxation reads constants too where it is formed at the samples
  localKernelPrototype_.bindGlobals(*global.onHost);

  fsgFlux_.bindGlobals(*global.onHost);
  dirichletFlux_.bindGlobals(*global.onHost);
  nodalLfKrnlPrototype_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceVolumeKernelPrototype_.bindGlobals(*global.onDevice);
  deviceLocalFluxKernelPrototype_.bindGlobals(*global.onDevice);
  deviceFluxLocalAllKernelPrototype_.bindGlobals(*global.onDevice);
  deviceLocalKernelPrototype_.bindGlobals(*global.onDevice);
#endif
}

void Local::computeIntegral(
    real* timeIntegratedDoFs, LTS::Ref& data, LocalTmp& tmp, double time, double timeStepWidth) {
  // assert alignments
#ifndef NDEBUG
  assert((reinterpret_cast<uintptr_t>(timeIntegratedDoFs)) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(tmp.timeIntegratedAne)) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(data.get<LTS::Dofs>())) % Vectorsize == 0);
#endif

  alignas(Alignment) real Qext[tensor::Qext::size()];

  kernel::volumeExt volKrnl = volumeKernelPrototype_;
  volKrnl.Qext = Qext;
  volKrnl.I = timeIntegratedDoFs;
  kernels::bindStarOperands(volKrnl, data.get<LTS::LocalIntegration>());

  kernel::localFluxExt lfKrnl = localFluxKernelPrototype_;
  lfKrnl.Qext = Qext;
  lfKrnl.I = timeIntegratedDoFs;
  lfKrnl._prefetch.I = timeIntegratedDoFs + tensor::I::size();
  lfKrnl._prefetch.Q = data.get<LTS::Dofs>() + tensor::Q::size();

  volKrnl.execute();

  const auto& cellBoundaryMapping = data.get<LTS::BoundaryMapping>();
  const auto& materialData = data.get<LTS::Material>();

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // no element local contribution in the case of dynamic rupture boundary conditions
    if (data.get<LTS::CellInformation>().faceTypes[face] != FaceType::DynamicRupture) {
      kernels::bindLocalFluxOperands(lfKrnl, data.get<LTS::LocalIntegration>(), face);
      lfKrnl.execute(face);
    }

    switch (data.get<LTS::CellInformation>().faceTypes[face]) {
    case FaceType::FreeSurfaceGravity: {
      auto kernel = fsgFlux_;
      kernel.g2m = -2 * this->gravitationalAcceleration_;

      const real localRho = materialData.local->getDensity();
      kernel.rho = &localRho;
      kernel.averageNormalDisplacement = tmp.nodalAvgDisplacements[face].data();

      kernel.Qext = Qext;
      kernel.AminusT = data.get<LTS::NeighboringIntegration>().nAmNm1[face];

      kernel.execute(face);
      break;
    }
    case FaceType::Dirichlet: {
      auto kernel = dirichletFlux_;
      kernel.dirichletOffset = cellBoundaryMapping[face].dirichletOffset;
      kernel.dt = timeStepWidth;

      kernel.Qext = Qext;
      kernel.AminusT = data.get<LTS::NeighboringIntegration>().nAmNm1[face];

      kernel.execute(face);
      break;
    }
    case FaceType::Analytical: {
      assert(initConds_ != nullptr);
      const auto applyAnalyticalSolution = kernels::ApplyAnalyticalSolution(initConds_, data);

      alignas(Alignment) real dofsFaceBoundaryNodal[tensor::INodal::size()];
      analyticalBoundary_.evaluate(cellBoundaryMapping[face],
                                   applyAnalyticalSolution,
                                   dofsFaceBoundaryNodal,
                                   time,
                                   timeStepWidth);

      auto nodalLfKrnl = nodalLfKrnlPrototype_;
      nodalLfKrnl.Qext = Qext;
      nodalLfKrnl.INodal = dofsFaceBoundaryNodal;
      nodalLfKrnl.AminusT = data.get<LTS::NeighboringIntegration>().nAmNm1[face];
      nodalLfKrnl.execute(face);
      break;
    }
    default:
      // No boundary condition.
      break;
    }
  }

  kernel::local lKrnl = localKernelPrototype_;
  // where the material varies inside the cell, the relaxation is formed
  // from what it says at the sample points and the kernel takes no
  // matrix at all
  set_E(lKrnl, data.get<LTS::LocalIntegration>().specific.E);
  kernels::bindSourceOperands(lKrnl, data.get<LTS::LocalIntegration>());
  lKrnl.Iane = tmp.timeIntegratedAne;
  lKrnl.Q = data.get<LTS::Dofs>();
  lKrnl.Qane = data.get<LTS::DofsAne>();
  lKrnl.Qext = Qext;
  lKrnl.W = data.get<LTS::LocalIntegration>().specific.W;
  lKrnl.w = data.get<LTS::LocalIntegration>().specific.w;

  lKrnl.execute();
}

PerformanceEstimate Local::metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes) const {
  auto estimate = PerformanceEstimate::fromKernel<seissol::kernel::volumeExt>();

#if defined(ACL_DEVICE) && defined(SEISSOL_DEVICE_COMBINE_LOCAL_FLUX)
  constexpr bool CombineLocalFlux = true;
#else
  constexpr bool CombineLocalFlux = false;
#endif

  if constexpr (CombineLocalFlux) {
    estimate += PerformanceEstimate::fromKernel<seissol::kernel::fluxLocalAll>();
  } else {
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (faceTypes[face] != FaceType::DynamicRupture) {
        estimate += PerformanceEstimate::fromKernel<seissol::kernel::localFluxExt>(face);
      }
      switch (faceTypes[face]) {
      case FaceType::FreeSurfaceGravity:
        estimate += PerformanceEstimate::fromKernel<seissol::kernel::fsgFlux>(face);
        break;
      case FaceType::Dirichlet:
        estimate += PerformanceEstimate::fromKernel<seissol::kernel::dirichletFlux>(face);
        break;
      case FaceType::Analytical:
        estimate += PerformanceEstimate::fromKernel<seissol::kernel::localFluxNodal>(face);
        break;
      default:
        break;
      }
    }

    estimate += PerformanceEstimate::fromKernel<seissol::kernel::local>();
  }

  // legacy memory estimate
  std::uint64_t reals = 0;

  // star matrices load
  reals += yateto::computeFamilySize<tensor::star>() + tensor::w::size() + tensor::W::size() +
           tensor::E::size();
  // flux solvers
  reals += 4 * tensor::AplusT::size();

  // DOFs write
  reals += tensor::Q::size() + tensor::Qane::size();

  estimate.bytes = reals * sizeof(real);

  return estimate;
}

void Local::computeBatchedIntegral(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& dataTable,
    SEISSOL_GPU_PARAM recording::ConditionalIndicesTable& indicesTable,
    SEISSOL_GPU_PARAM double timeStepWidth,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;
  // Volume integral
  ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  kernel::gpu_volumeExt volKrnl = deviceVolumeKernelPrototype_;

  if (dataTable.find(key) != dataTable.end()) {
    auto& entry = dataTable[key];

    volKrnl.numElements = (dataTable[key].get(inner_keys::Wp::Id::Dofs))->getSize();

    volKrnl.I =
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
    volKrnl.Qext = (entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();

    kernels::bindStarOperandsBatched(
        volKrnl,
        const_cast<const real**>(
            (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr()));
    volKrnl.streamPtr = runtime.stream();
    volKrnl.execute();

#ifdef SEISSOL_DEVICE_COMBINE_LOCAL_FLUX
    auto krnl = deviceFluxLocalAllKernelPrototype_;

    krnl.numElements = entry.get(inner_keys::Wp::Id::Dofs)->getSize();
    krnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    krnl.Qane = (entry.get(inner_keys::Wp::Id::DofsAne))->getDeviceDataPtr();
    krnl.Qext = (entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();
    krnl.Iane =
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::IdofsAne))->getDeviceDataPtr());
    krnl.W = const_cast<const real**>(
        entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    krnl.extraOffset_W = SEISSOL_OFFSET(LocalIntegrationData, specific.W);
    krnl.w = const_cast<const real**>(
        entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    krnl.extraOffset_w = SEISSOL_OFFSET(LocalIntegrationData, specific.w);
    // where the material varies inside the cell, the relaxation is formed
    // from what it says at the sample points and the kernel takes no
    // matrix at all
    set_E(krnl,
          const_cast<const real**>(
              entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr()));
    set_extraOffset_E(krnl, SEISSOL_OFFSET(LocalIntegrationData, specific.E));
    kernels::bindSourceOperandsBatched(
        krnl,
        const_cast<const real**>(
            entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr()));
    krnl.streamPtr = runtime.stream();

    SEISSOL_OFFSET_ASSERT(LocalIntegrationData, specific.W);
    SEISSOL_OFFSET_ASSERT(LocalIntegrationData, specific.w);
    SEISSOL_OFFSET_ASSERT(LocalIntegrationData, specific.E);

    krnl.I = const_cast<const real**>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());

    kernels::bindLocalFluxAllOperandsBatched(
        krnl,
        const_cast<const real**>(
            entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr()));

    krnl.execute();
#endif
  }

// deprecated code; kept for comparison reasons
#ifndef SEISSOL_DEVICE_COMBINE_LOCAL_FLUX
  kernel::gpu_localFluxExt localFluxKrnl = deviceLocalFluxKernelPrototype_;
  kernel::gpu_local localKrnl = deviceLocalKernelPrototype_;

  // Local Flux Integral
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    key = ConditionalKey(*KernelNames::LocalFlux, !FaceKinds::DynamicRupture, face);

    if (dataTable.find(key) != dataTable.end()) {
      auto& entry = dataTable[key];
      localFluxKrnl.numElements = entry.get(inner_keys::Wp::Id::Dofs)->getSize();
      localFluxKrnl.Qext = (entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();
      localFluxKrnl.I =
          const_cast<const real**>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
      kernels::bindLocalFluxOperandsBatched(
          localFluxKrnl,
          const_cast<const real**>(
              entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr()),
          face);
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
    set_E(localKrnl,
          const_cast<const real**>(
              entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr()));
    set_extraOffset_E(localKrnl, SEISSOL_OFFSET(LocalIntegrationData, specific.E));
    kernels::bindSourceOperandsBatched(
        localKrnl,
        const_cast<const real**>(
            entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr()));
    localKrnl.streamPtr = runtime.stream();

    SEISSOL_OFFSET_ASSERT(LocalIntegrationData, specific.W);
    SEISSOL_OFFSET_ASSERT(LocalIntegrationData, specific.w);
    SEISSOL_OFFSET_ASSERT(LocalIntegrationData, specific.E);

    localKrnl.execute();
  }
#endif

#else
  logError() << "No GPU implementation provided";
#endif
}

void Local::evaluateBatchedTimeDependentBc(recording::ConditionalPointersToRealsTable& dataTable,
                                           recording::ConditionalIndicesTable& indicesTable,
                                           LTS::Layer& layer,
                                           double time,
                                           double timeStepWidth,
                                           seissol::parallel::runtime::StreamRuntime& runtime) {}

} // namespace seissol::kernels::solver::linearckanelastic
