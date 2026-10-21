// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "Kernels/LinearCK/Local.h"

#include "Alignment.h"
#include "Common/Constants.h"
#include "Common/Marker.h"
#include "Common/Typedefs.h"
#include "Config.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/Typedefs.h"
#include "Kernels/AnalyticalBoundary.h"
#include "Kernels/Common.h"
#include "Kernels/Interface.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/LTS.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdint.h>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"

#include <utility>
#endif

#ifndef ACL_DEVICE
#include <utils/logger.h>
#endif

GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(extraOffset_ET)
GENERATE_HAS_MEMBER(sourceMatrix)

namespace seissol::kernels::solver::linearck {

template <typename Cfg>
void Local<Cfg>::setGlobalData(const CompoundGlobalData<Cfg>& global) {
  volumeKernelPrototype_.bindGlobals(*global.onHost);
  localFluxKernelPrototype_.bindGlobals(*global.onHost);
  nodalLfKrnlPrototype_.bindGlobals(*global.onHost);
  ghostStateBoundary_.bindGlobals(*global.onHost);
  fsgFlux_.bindGlobals(*global.onHost);
  dirichletFlux_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceVolumeKernelPrototype_.bindGlobals(*global.onDevice);
  deviceLocalFluxKernelPrototype_.bindGlobals(*global.onDevice);
  deviceLocalFluxAllKernelPrototype_.bindGlobals(*global.onDevice);
  deviceNodalLfKrnlPrototype_.bindGlobals(*global.onDevice);
  deviceFsgFlux_.bindGlobals(*global.onDevice);
  deviceDirichletFlux_.bindGlobals(*global.onDevice);
  deviceGhostStates_ = std::make_shared<DeviceGhostStates<Cfg>>(global.onDevice);
#endif
}

template <typename Cfg>
void Local<Cfg>::computeIntegral(real* timeIntegratedDoFs,
                                 LTS::Ref<Cfg>& data,
                                 LocalTmp<Cfg>& tmp,
                                 double time,
                                 double timeStepWidth) {
  assert(reinterpret_cast<uintptr_t>(timeIntegratedDoFs) % Vectorsize == 0);
  assert(reinterpret_cast<uintptr_t>(data.template get<LTS::Dofs>()) % Vectorsize == 0);

  const auto& materialData = data.template get<LTS::Material>();
  const auto& cellBoundaryMapping = data.template get<LTS::BoundaryMapping>();

  // the kernels add their products into Q term by term; into the DOFs, every one of these
  // additions would round at the size of the DOFs rather than of the update, so that the rounding
  // error of a step grows with the number of terms. So they go into a buffer, which is added once.
  alignas(Alignment) real update[tensor::Q<Cfg>::size()]{};

  kernel::volume<Cfg> volKrnl = volumeKernelPrototype_;
  volKrnl.Q = update;
  volKrnl.I = timeIntegratedDoFs;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
    volKrnl.star(i) = data.template get<LTS::LocalIntegration>().starMatrices[i];
  }

  // Optional source term
  set_ET(volKrnl, get_ptr_sourceMatrix(data.template get<LTS::LocalIntegration>().specific));

  kernel::localFlux<Cfg> lfKrnl = localFluxKernelPrototype_;
  lfKrnl.Q = update;
  lfKrnl.I = timeIntegratedDoFs;
  lfKrnl._prefetch.I = timeIntegratedDoFs + tensor::I<Cfg>::size();
  lfKrnl._prefetch.Q = data.template get<LTS::Dofs>() + tensor::Q<Cfg>::size();

  volKrnl.execute();

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // no element local contribution in the case of dynamic rupture boundary conditions
    if (data.template get<LTS::CellInformation>().faceTypes[face] != FaceType::DynamicRupture) {
      lfKrnl.AplusT = data.template get<LTS::LocalIntegration>().nApNm1[face];
      lfKrnl.execute(face);
    }

    alignas(Alignment) real dofsFaceBoundaryNodal[tensor::INodal<Cfg>::size()];
    auto nodalLfKrnl = nodalLfKrnlPrototype_;
    nodalLfKrnl.Q = update;
    nodalLfKrnl.INodal = dofsFaceBoundaryNodal;
    nodalLfKrnl._prefetch.I = timeIntegratedDoFs + tensor::I<Cfg>::size();
    nodalLfKrnl._prefetch.Q = data.template get<LTS::Dofs>() + tensor::Q<Cfg>::size();
    nodalLfKrnl.AminusT = data.template get<LTS::NeighboringIntegration>().nAmNm1[face];

    // Include some boundary conditions here.
    switch (data.template get<LTS::CellInformation>().faceTypes[face]) {
    case FaceType::FreeSurfaceGravity: {
      auto kernel = fsgFlux_;
      kernel.g2m = -2 * this->gravitationalAcceleration_;

      const real localRho = materialData.local->getDensity();
      kernel.rho = &localRho;
      kernel.averageNormalDisplacement = tmp.nodalAvgDisplacements[face].data();

      kernel.Q = update;
      kernel.AminusT = data.template get<LTS::NeighboringIntegration>().nAmNm1[face];

      kernel.execute(face);
      break;
    }
    case FaceType::Dirichlet: {
      auto* dirichletOffset = cellBoundaryMapping[face].dirichletOffset;
      assert(dirichletOffset != nullptr);

      auto kernel = dirichletFlux_;
      kernel.dirichletOffset = dirichletOffset;
      kernel.dt = timeStepWidth;

      kernel.Q = update;
      kernel.AminusT = data.template get<LTS::NeighboringIntegration>().nAmNm1[face];

      kernel.execute(face);
      break;
    }
    case FaceType::Analytical: {
      ghostStateBoundary_.analytical(this->analyticalGhostState(),
                                     this->initConds_,
                                     data,
                                     face,
                                     time,
                                     timeStepWidth,
                                     dofsFaceBoundaryNodal);
      nodalLfKrnl.execute(face);
      break;
    }
    case FaceType::NonlinearDirichlet: {
      assert(this->nonlinearDirichlet() != nullptr);
      assert(tmp.timeDerivatives != nullptr);
      // the state of the cell at a time into the step: from its Taylor series, or from the
      // polynomial in time of the space-time predictor
      const real* derivatives = tmp.timeDerivatives;
      const auto stateAt = [derivatives, timeStepWidth](double tau, real* state) {
        // the predictor reads it, the Taylor series does not
        static_cast<void>(timeStepWidth);
        if constexpr (Cfg::Solver == SolverType::LinearCK) {
          kernel::derivativeTaylorExpansion<Cfg> taylorKrnl;
          taylorKrnl.I = state;
          double power = 1;
          for (std::size_t d = 0; d < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++d) {
            taylorKrnl.dQ(d) = derivatives + yateto::computeFamilySize<tensor::dQ<Cfg>>(1, d);
            taylorKrnl.power(d) = static_cast<real>(power);
            power *= tau / static_cast<double>(d + 1);
          }
          taylorKrnl.execute();
        } else {
#ifdef SEISSOL_KERNELS_STP
          const auto basis = seissol::kernels::timeBasis<Cfg>().point(tau, timeStepWidth);
          kernel::evaluateDOFSAtTimeSTP<Cfg> evaluateKrnl;
          evaluateKrnl.spaceTimePredictor = derivatives;
          evaluateKrnl.timeBasisFunctionsAtPoint = basis.data();
          evaluateKrnl.QAtTimeSTP = state;
          evaluateKrnl.execute();
#endif
        }
      };
      ghostStateBoundary_.nonlinearDirichlet(*this->nonlinearDirichlet(),
                                             stateAt,
                                             face,
                                             cellBoundaryMapping[face],
                                             data.template get<LTS::MaterialParameters>(),
                                             time,
                                             timeStepWidth,
                                             dofsFaceBoundaryNodal);
      nodalLfKrnl.execute(face);
      break;
    }
    default:
      // No boundary condition.
      break;
    }
  }

  real* dofs = data.template get<LTS::Dofs>();
#pragma omp simd
  for (std::size_t i = 0; i < tensor::Q<Cfg>::size(); ++i) {
    dofs[i] += update[i];
  }
}

template <typename Cfg>
void Local<Cfg>::computeBatchedIntegral(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& dataTable,
    SEISSOL_GPU_PARAM recording::ConditionalIndicesTable& /*indicesTable*/,
    SEISSOL_GPU_PARAM double timeStepWidth,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;
  // Volume integral
  const ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  kernel::gpu_volume<Cfg> volKrnl = deviceVolumeKernelPrototype_;
  const kernel::gpu_localFlux<Cfg> localFluxKrnl = deviceLocalFluxKernelPrototype_;

  // the boundary kernels below run on the same temporary memory
  const auto maxTmpMem =
      yateto::getMaxTmpMemRequired(volKrnl, localFluxKrnl, deviceFsgFlux_, deviceDirichletFlux_);

  // volume kernel always contains more elements than any local one
  const auto maxNumElements = dataTable.find(key) != dataTable.end()
                                  ? (dataTable[key].get<real*>(inner_keys::Wp::Id::Dofs))->getSize()
                                  : 0;
  auto tmpMem = runtime.memoryHandle<real>((maxTmpMem * maxNumElements) / sizeof(real));
  if (dataTable.find(key) != dataTable.end()) {
    auto& entry = dataTable[key];

    volKrnl.numElements = maxNumElements;

    volKrnl.Q = (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    volKrnl.I =
        const_cast<const real**>((entry.get<real*>(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());

    const auto** localIntegrationPtrs = const_cast<const real**>(
        (entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());

    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData<Cfg>, starMatrices);
    for (size_t i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
      volKrnl.star(i) = localIntegrationPtrs;
      volKrnl.extraOffset_star(i) =
          SEISSOL_ARRAY_OFFSET(LocalIntegrationData<Cfg>, starMatrices, i);
    }

    constexpr auto SourceMatrixOffset =
        offsetof(LocalIntegrationData<Cfg>, specific) +
        get_offset_sourceMatrix<decltype(LocalIntegrationData<Cfg>::specific)>();
    static_assert(SourceMatrixOffset % sizeof(real) == 0,
                  "SourceMatrixOffset is not dividable by the real size.");

    set_ET(volKrnl, localIntegrationPtrs);
    set_extraOffset_ET(volKrnl, SourceMatrixOffset / sizeof(real));

    volKrnl.linearAllocator.initialize(tmpMem.get());
    volKrnl.streamPtr = runtime.stream();
    volKrnl.execute();

#ifdef SEISSOL_DEVICE_COMBINE_LOCAL_FLUX
    kernel::gpu_localFluxAll<Cfg> localFluxKrnl = deviceLocalFluxAllKernelPrototype_;
    localFluxKrnl.numElements = entry.get<real*>(inner_keys::Wp::Id::Dofs)->getSize();
    localFluxKrnl.Q = (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    localFluxKrnl.I =
        const_cast<const real**>((entry.get<real*>(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());

    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData<Cfg>, nApNm1);
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      localFluxKrnl.AplusTAll(face) = const_cast<const real**>(
          entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
      localFluxKrnl.extraOffset_AplusTAll(face) =
          SEISSOL_ARRAY_OFFSET(LocalIntegrationData<Cfg>, nApNm1, face);
    }
    localFluxKrnl.linearAllocator.initialize(tmpMem.get());
    localFluxKrnl.streamPtr = runtime.stream();
    localFluxKrnl.execute();
#endif
  }

  // Local Flux Integral
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // deprecated code; kept for comparison reasons
#ifndef SEISSOL_DEVICE_COMBINE_LOCAL_FLUX
    const auto key = ConditionalKey(*KernelNames::LocalFlux, !FaceKinds::DynamicRupture, face);

    if (dataTable.find(key) != dataTable.end()) {
      auto& entry = dataTable[key];
      localFluxKrnl.numElements = entry.get<real*>(inner_keys::Wp::Id::Dofs)->getSize();
      localFluxKrnl.Q = (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
      localFluxKrnl.I = const_cast<const real**>(
          (entry.get<real*>(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
      localFluxKrnl.AplusT = const_cast<const real**>(
          entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());

      SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData<Cfg>, nApNm1);
      localFluxKrnl.extraOffset_AplusT =
          SEISSOL_ARRAY_OFFSET(LocalIntegrationData<Cfg>, nApNm1, face);
      localFluxKrnl.linearAllocator.initialize(tmpMem.get());
      localFluxKrnl.streamPtr = runtime.stream();
      localFluxKrnl.execute(face);
    }
#endif

    const ConditionalKey fsgKey(
        *KernelNames::BoundaryConditions, *ComputationKind::FreeSurfaceGravity, face);
    if (dataTable.find(fsgKey) != dataTable.end()) {
      auto** nodalAvgDisplacements = dataTable[fsgKey]
                                         .get<real*>(inner_keys::Wp::Id::NodalAvgDisplacements)
                                         ->getDeviceDataPtr();
      auto** rhos = dataTable[fsgKey].get<real*>(inner_keys::Wp::Id::FSGData)->getDeviceDataPtr();

      auto bcKernel = deviceFsgFlux_;
      bcKernel.g2m = -2 * this->gravitationalAcceleration_;
      bcKernel.rho = const_cast<const real**>(rhos);
      bcKernel.extraOffset_rho = 2;
      bcKernel.averageNormalDisplacement = const_cast<const real**>(nodalAvgDisplacements);
      bcKernel.Q = (dataTable[fsgKey].get<real*>(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
      bcKernel.AminusT =
          const_cast<const real**>(dataTable[fsgKey]
                                       .get<real*>(inner_keys::Wp::Id::NeighborIntegrationData)
                                       ->getDeviceDataPtr());
      bcKernel.extraOffset_AminusT =
          SEISSOL_ARRAY_OFFSET(NeighboringIntegrationData<Cfg>, nAmNm1, face);

      bcKernel.numElements = dataTable[fsgKey].get<real*>(inner_keys::Wp::Id::Dofs)->getSize();

      bcKernel.linearAllocator.initialize(tmpMem.get());
      bcKernel.streamPtr = runtime.stream();

      bcKernel.execute(face);
    }

    const ConditionalKey dirichletKey(
        *KernelNames::BoundaryConditions, *ComputationKind::Dirichlet, face);
    if (dataTable.find(dirichletKey) != dataTable.end()) {
      auto* dirichletOffsetPtrs = dataTable[dirichletKey]
                                      .get<real*>(inner_keys::Wp::Id::DirichletOffset)
                                      ->getDeviceDataPtr();

      auto bcKernel = deviceDirichletFlux_;
      bcKernel.dirichletOffset = const_cast<const real**>(dirichletOffsetPtrs);
      bcKernel.dt = timeStepWidth;
      bcKernel.Q =
          (dataTable[dirichletKey].get<real*>(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
      bcKernel.AminusT =
          const_cast<const real**>(dataTable[dirichletKey]
                                       .get<real*>(inner_keys::Wp::Id::NeighborIntegrationData)
                                       ->getDeviceDataPtr());
      bcKernel.extraOffset_AminusT =
          SEISSOL_ARRAY_OFFSET(NeighboringIntegrationData<Cfg>, nAmNm1, face);

      bcKernel.numElements =
          dataTable[dirichletKey].get<real*>(inner_keys::Wp::Id::Dofs)->getSize();

      bcKernel.linearAllocator.initialize(tmpMem.get());
      bcKernel.streamPtr = runtime.stream();

      bcKernel.execute(face);
    }
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

template <typename Cfg>
void Local<Cfg>::evaluateBatchedTimeDependentBc(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& dataTable,
    SEISSOL_GPU_PARAM recording::ConditionalIndicesTable& indicesTable,
    SEISSOL_GPU_PARAM LTS::Layer& layer,
    SEISSOL_GPU_PARAM double time,
    SEISSOL_GPU_PARAM double timeStepWidth,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {

#ifdef ACL_DEVICE
  using namespace seissol::recording;

  // the flux of the ghost states of the faces of a batch, from their slots of AnalyticScratch
  const auto nodalFlux = [&](const ConditionalKey& key, std::size_t face) {
    const std::size_t numElements = dataTable[key].get<real*>(inner_keys::Wp::Id::Dofs)->getSize();
    auto nodalLfKrnl = deviceNodalLfKrnlPrototype_;
    auto tmpMem = runtime.memoryHandle<real>(
        (yateto::getMaxTmpMemRequired(nodalLfKrnl) * numElements) / sizeof(real));
    nodalLfKrnl.linearAllocator.initialize(tmpMem.get());
    nodalLfKrnl.INodal = const_cast<const real**>(
        dataTable[key].get<real*>(inner_keys::Wp::Id::Analytical)->getDeviceDataPtr());
    nodalLfKrnl.AminusT = const_cast<const real**>(
        dataTable[key].get<real*>(inner_keys::Wp::Id::NeighborIntegrationData)->getDeviceDataPtr());
    nodalLfKrnl.extraOffset_AminusT =
        SEISSOL_ARRAY_OFFSET(NeighboringIntegrationData<Cfg>, nAmNm1, face);
    nodalLfKrnl.Q = dataTable[key].get<real*>(inner_keys::Wp::Id::Dofs)->getDeviceDataPtr();
    nodalLfKrnl.streamPtr = runtime.stream();
    nodalLfKrnl.numElements = numElements;
    nodalLfKrnl.execute(face);
  };

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    const ConditionalKey analyticalKey(
        *KernelNames::BoundaryConditions, *ComputationKind::Analytical, face);
    if (indicesTable.find(analyticalKey) != indicesTable.end()) {
      // the setup refuses an initial condition that is no script that compiles in a device build
      assert(this->analyticalGhostState() != nullptr);
      deviceGhostStates_->analytical(*this->analyticalGhostState(),
                                     indicesTable,
                                     layer,
                                     face,
                                     time,
                                     timeStepWidth,
                                     runtime.stream());
      nodalFlux(analyticalKey, face);
    }

    const ConditionalKey nonlinearDirichletKey(
        *KernelNames::BoundaryConditions, *ComputationKind::NonlinearDirichlet, face);
    if (indicesTable.find(nonlinearDirichletKey) != indicesTable.end()) {
      assert(this->nonlinearDirichlet() != nullptr);
      deviceGhostStates_->nonlinearDirichlet(*this->nonlinearDirichlet(),
                                             dataTable,
                                             indicesTable,
                                             layer,
                                             face,
                                             time,
                                             timeStepWidth,
                                             runtime);
      nodalFlux(nonlinearDirichletKey, face);
    }
  }
#else
  logError() << "No GPU implementation provided";
#endif // ACL_DEVICE
}

template <typename Cfg>
PerformanceEstimate
    Local<Cfg>::metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes) const {
  PerformanceEstimate estimate;
  estimate += PerformanceEstimate::fromKernel<seissol::kernel::volume<Cfg>>();

#if defined(ACL_DEVICE) && defined(SEISSOL_DEVICE_COMBINE_LOCAL_FLUX)
  constexpr bool CombineLocalFlux = true;
#else
  constexpr bool CombineLocalFlux = false;
#endif

  if constexpr (CombineLocalFlux) {
    estimate += PerformanceEstimate::fromKernel<seissol::kernel::localFluxAll<Cfg>>();
  }

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // Local flux is executed for all faces that are not dynamic rupture.
    // For those cells, the flux is taken into account during the neighbor kernel.
    if (faceTypes[face] != FaceType::DynamicRupture && !CombineLocalFlux) {
      estimate += PerformanceEstimate::fromKernel<seissol::kernel::localFlux<Cfg>>(face);
    }

    // Take boundary condition flops into account.
    // Note that this only includes the flops of the kernels but not of the
    // boundary condition implementation.
    // The (probably incorrect) assumption is that they are negligible.
    switch (faceTypes[face]) {
    case FaceType::FreeSurfaceGravity:
      estimate += PerformanceEstimate::fromKernel<seissol::kernel::fsgFlux<Cfg>>(face);
      break;
    case FaceType::Dirichlet:
      estimate += PerformanceEstimate::fromKernel<seissol::kernel::dirichletFlux<Cfg>>(face);
      break;
    case FaceType::Analytical:
      // and the script, at every node, integrated over the step
      estimate += PerformanceEstimate::fromKernel<seissol::kernel::localFluxNodal<Cfg>>(face);
      break;
    case FaceType::NonlinearDirichlet:
      // and the script, at every node, integrated over the step
      estimate += PerformanceEstimate::fromKernel<seissol::kernel::localFluxNodal<Cfg>>(face);
      if constexpr (Cfg::Solver == SolverType::LinearCK) {
        estimate +=
            (PerformanceEstimate::fromKernel<seissol::kernel::projectToFaceNodes<Cfg>>(face) +
             PerformanceEstimate::fromKernel<seissol::kernel::derivativeTaylorExpansion<Cfg>>()) *
            Cfg::ConvergenceOrder;
      } else {
#ifdef SEISSOL_KERNELS_STP
        estimate +=
            (PerformanceEstimate::fromKernel<seissol::kernel::projectToFaceNodes<Cfg>>(face) +
             PerformanceEstimate::fromKernel<seissol::kernel::evaluateDOFSAtTimeSTP<Cfg>>()) *
            Cfg::ConvergenceOrder;
#endif
      }
      break;
    default:
      break;
    }
  }

  // legacy memory estimate
  std::uint64_t reals = 0;

  // star matrices load
  reals += yateto::computeFamilySize<tensor::star<Cfg>>();
  // flux solvers
  reals += static_cast<std::uint64_t>(4 * tensor::AplusT<Cfg>::size());

  // DOFs write
  reals += tensor::Q<Cfg>::size();

  estimate.bytes = reals * sizeof(real);

  return estimate;
}

#define SEISSOL_INSTANTIATE(Cfg) template class Local<Cfg>;
SEISSOL_FOR_EACH_CONFIG_LINEARCK(SEISSOL_INSTANTIATE)
SEISSOL_FOR_EACH_CONFIG_STP(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels::solver::linearck
