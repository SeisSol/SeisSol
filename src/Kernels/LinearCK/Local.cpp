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
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/BoundarySetup.h"
#include "Initializer/Typedefs.h"
#include "Kernels/AnalyticalBoundary.h"
#include "Kernels/Common.h"
#include "Kernels/FaceTraces.h"
#include "Kernels/Interface.h"
#include "Kernels/Precision.h"
#include "Kernels/StarOperands.h"
#include "Kernels/WeightAdjustment.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/OperatorLayout.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdint.h>
#include <utils/logger.h>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#endif

GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)
GENERATE_HAS_MEMBER(volumeFused)

namespace seissol::kernels::solver::linearck {

static_assert(!WeightAdjustable || kernels::size<tensor::fluxBuffer>() == tensor::Q::size(),
              "the faces of a weight-adjusted cell lift into its flux buffer as into its field");

namespace {
/// Hands the volume kernel what it reads: the volume term the predictor of the cell formed
/// already, where the build forms it there (LocalTmp::volumeFused), or else the time integral
/// and the operator of the cell.
template <typename KernelT>
void bindVolumeOperands(KernelT& krnl, real* timeIntegratedDoFs, LTS::Ref& data, LocalTmp& tmp) {
  if constexpr (has_volumeFused<KernelT>::value) {
    krnl.volumeFused = tmp.volumeFused;
  } else {
    krnl.I = timeIntegratedDoFs;
    kernels::bindStarOperands(krnl, data.get<LTS::LocalIntegration>());

    // Optional source term
    set_ET(krnl, get_ptr_sourceMatrix(data.get<LTS::LocalIntegration>().specific));
    kernels::bindSourceOperands(krnl, data.get<LTS::LocalIntegration>());
  }
}
} // namespace

void Local::setGlobalData(const CompoundGlobalData& global) {
  volumeKernelPrototype_.bindGlobals(*global.onHost);
  faceTraceKernelPrototype_.bindGlobals(*global.onHost);
  nodalLfKrnlPrototype_.bindGlobals(*global.onHost);
  fsgFlux_.bindGlobals(*global.onHost);
  dirichletFlux_.bindGlobals(*global.onHost);
  foldedBoundaryFlux_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceFaceTracesKernelPrototype_.bindGlobals(*global.onDevice);
  deviceNodalLfKrnlPrototype_.bindGlobals(*global.onDevice);
  deviceFsgFlux_.bindGlobals(*global.onDevice);
  deviceDirichletFlux_.bindGlobals(*global.onDevice);
  deviceFoldedBoundaryFlux_.bindGlobals(*global.onDevice);
#endif
}

void Local::computeIntegral(
    real* timeIntegratedDoFs, LTS::Ref& data, LocalTmp& tmp, double time, double timeStepWidth) {
  assert(reinterpret_cast<uintptr_t>(timeIntegratedDoFs) % Vectorsize == 0);
  assert(reinterpret_cast<uintptr_t>(data.get<LTS::Dofs>()) % Vectorsize == 0);

  const auto& materialData = data.get<LTS::Material>();
  const auto& cellBoundaryMapping = data.get<LTS::BoundaryMapping>();

  // A weight-adjusted cell lifts its face terms without material into its flux buffer, which
  // starts empty here and which the neighbour integral finishes (WeightAdjustable); the volume
  // term goes to the field either way.
  real* fluxTarget = data.get<LTS::Dofs>();
  if (WeightAdjustable && data.get<LTS::CellInformation>().weightAdjusted) {
    fluxTarget = data.get<LTS::FluxBuffer>();
    std::fill_n(fluxTarget, kernels::size<tensor::fluxBuffer>(), static_cast<real>(0));
  }

  kernel::volume volKrnl = volumeKernelPrototype_;
  volKrnl.Q = data.get<LTS::Dofs>();
  bindVolumeOperands(volKrnl, timeIntegratedDoFs, data, tmp);

  volKrnl.execute();

  // the trace of the time integral on every face, which the corrector of the cell reads together
  // with the one its neighbor writes there
  kernels::writeFaceTraces(faceTraceKernelPrototype_, timeIntegratedDoFs, data);

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    alignas(Alignment) real dofsFaceBoundaryNodal[tensor::INodal::size()];
    auto nodalLfKrnl = nodalLfKrnlPrototype_;
    nodalLfKrnl.Q = fluxTarget;
    nodalLfKrnl.INodal = dofsFaceBoundaryNodal;
    nodalLfKrnl._prefetch.I = timeIntegratedDoFs + tensor::I::size();
    nodalLfKrnl._prefetch.Q = data.get<LTS::Dofs>() + tensor::Q::size();
    nodalLfKrnl.AminusT = data.get<LTS::NeighboringIntegration>().nAmNm1[face];

    // what the condition of the face folds into the local flux solver, where the corrector does
    // not apply it
    if (appliesFoldApart(data.get<LTS::CellInformation>().faceTypes[face])) {
      auto kernel = foldedBoundaryFlux_;
      kernel.Q = fluxTarget;
      kernel.I = timeIntegratedDoFs;
      kernel.AplusT = data.get<LTS::LocalIntegration>().nApNm1[face];
      kernel.execute(face);
    }

    // Include some boundary conditions here.
    switch (data.get<LTS::CellInformation>().faceTypes[face]) {
    case FaceType::FreeSurfaceGravity: {
      auto kernel = fsgFlux_;
      kernel.g2m = -2 * this->gravitationalAcceleration_;

      const real localRho = materialData.local->getDensity();
      kernel.rho = &localRho;
      kernel.averageNormalDisplacement = tmp.nodalAvgDisplacements[face].data();

      kernel.Q = fluxTarget;
      kernel.AminusT = data.get<LTS::NeighboringIntegration>().nAmNm1[face];

      kernel.execute(face);
      break;
    }
    case FaceType::Dirichlet: {
      auto* dirichletOffset = cellBoundaryMapping[face].dirichletOffset;
      assert(dirichletOffset != nullptr);

      auto kernel = dirichletFlux_;
      kernel.dirichletOffset = dirichletOffset;
      kernel.dt = timeStepWidth;

      kernel.Q = fluxTarget;
      kernel.AminusT = data.get<LTS::NeighboringIntegration>().nAmNm1[face];

      kernel.execute(face);
      break;
    }
    case FaceType::Analytical: {
      assert(initConds_ != nullptr);
      const auto applyAnalyticalSolution = ApplyAnalyticalSolution(initConds_, data);

      analyticalBoundary_.evaluate(cellBoundaryMapping[face],
                                   applyAnalyticalSolution,
                                   dofsFaceBoundaryNodal,
                                   time,
                                   timeStepWidth);
      nodalLfKrnl.execute(face);
      break;
    }
    default:
      // No boundary condition.
      break;
    }
  }
}

void Local::computeBatchedIntegral(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& dataTable,
    SEISSOL_GPU_PARAM recording::ConditionalIndicesTable& /*indicesTable*/,
    SEISSOL_GPU_PARAM double timeStepWidth,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;
  // The predictor has added the volume term to the field already (Spacetime::computeBatchedAder).
  const ConditionalKey key(KernelNames::Time || KernelNames::Volume);

  // the boundary kernels below run on the same temporary memory
  const auto maxTmpMem =
      yateto::getMaxTmpMemRequired(deviceFsgFlux_, deviceDirichletFlux_, deviceFoldedBoundaryFlux_);

  // the batch of the predictor has every cell, so at least as many as any batch of faces
  const auto maxNumElements = dataTable.find(key) != dataTable.end()
                                  ? (dataTable[key].get(inner_keys::Wp::Id::Dofs))->getSize()
                                  : 0;
  auto tmpMem = runtime.memoryHandle<real>((maxTmpMem * maxNumElements) / sizeof(real));

  // a weight-adjusted cell lifts its face terms into its flux buffer, which starts empty here and
  // which the neighbour integral finishes (WeightAdjustable); the flux kernels below write where
  // the batch says the faces of a cell lift into (inner_keys::Wp::Id::FluxTarget)
  WeightAdjustmentT::clearBatched(dataTable, runtime);

  // the trace of the time integral on every face, which the corrector of the cell reads together
  // with the one its neighbor writes there
  kernels::writeFaceTracesBatched(deviceFaceTracesKernelPrototype_, dataTable, runtime);

  // what the condition of the faces of a batch folds into the local flux solver, where the
  // corrector does not apply it (appliesFoldApart)
  const auto applyFoldApart = [&](auto& entry, std::size_t face) {
    auto foldKernel = deviceFoldedBoundaryFlux_;
    foldKernel.numElements = entry.get(inner_keys::Wp::Id::Dofs)->getSize();
    foldKernel.Q = entry.get(inner_keys::Wp::Id::FluxTarget)->getDeviceDataPtr();
    foldKernel.I =
        const_cast<const real**>(entry.get(inner_keys::Wp::Id::Idofs)->getDeviceDataPtr());
    foldKernel.AplusT = const_cast<const real**>(
        entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    foldKernel.extraOffset_AplusT = SEISSOL_ARRAY_OFFSET(LocalIntegrationData, nApNm1, face);
    foldKernel.linearAllocator.initialize(tmpMem.get());
    foldKernel.streamPtr = runtime.stream();
    foldKernel.execute(face);
  };

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    ConditionalKey fsgKey(
        *KernelNames::BoundaryConditions, *ComputationKind::FreeSurfaceGravity, face);
    if (dataTable.find(fsgKey) != dataTable.end()) {
      auto** nodalAvgDisplacements =
          dataTable[fsgKey].get(inner_keys::Wp::Id::NodalAvgDisplacements)->getDeviceDataPtr();
      auto** rhos = dataTable[fsgKey].get(inner_keys::Wp::Id::FSGData)->getDeviceDataPtr();

      auto bcKernel = deviceFsgFlux_;
      bcKernel.g2m = -2 * gravitationalAcceleration_;
      bcKernel.rho = const_cast<const real**>(rhos);
      bcKernel.extraOffset_rho = 2;
      bcKernel.averageNormalDisplacement = const_cast<const real**>(nodalAvgDisplacements);
      bcKernel.Q = (dataTable[fsgKey].get(inner_keys::Wp::Id::FluxTarget))->getDeviceDataPtr();
      bcKernel.AminusT = const_cast<const real**>(
          dataTable[fsgKey].get(inner_keys::Wp::Id::NeighborIntegrationData)->getDeviceDataPtr());
      bcKernel.extraOffset_AminusT = SEISSOL_ARRAY_OFFSET(NeighboringIntegrationData, nAmNm1, face);

      bcKernel.numElements = dataTable[fsgKey].get(inner_keys::Wp::Id::Dofs)->getSize();

      bcKernel.linearAllocator.initialize(tmpMem.get());
      bcKernel.streamPtr = runtime.stream();

      bcKernel.execute(face);

      if constexpr (appliesFoldApart(FaceType::FreeSurfaceGravity)) {
        applyFoldApart(dataTable[fsgKey], face);
      }
    }

    ConditionalKey dirichletKey(
        *KernelNames::BoundaryConditions, *ComputationKind::Dirichlet, face);
    if (dataTable.find(dirichletKey) != dataTable.end()) {
      auto* dirichletOffsetPtrs =
          dataTable[dirichletKey].get(inner_keys::Wp::Id::DirichletOffset)->getDeviceDataPtr();

      auto bcKernel = deviceDirichletFlux_;
      bcKernel.dirichletOffset = const_cast<const real**>(dirichletOffsetPtrs);
      bcKernel.dt = timeStepWidth;
      bcKernel.Q =
          (dataTable[dirichletKey].get(inner_keys::Wp::Id::FluxTarget))->getDeviceDataPtr();
      bcKernel.AminusT =
          const_cast<const real**>(dataTable[dirichletKey]
                                       .get(inner_keys::Wp::Id::NeighborIntegrationData)
                                       ->getDeviceDataPtr());
      bcKernel.extraOffset_AminusT = SEISSOL_ARRAY_OFFSET(NeighboringIntegrationData, nAmNm1, face);

      bcKernel.numElements = dataTable[dirichletKey].get(inner_keys::Wp::Id::Dofs)->getSize();

      bcKernel.linearAllocator.initialize(tmpMem.get());
      bcKernel.streamPtr = runtime.stream();

      bcKernel.execute(face);

      if constexpr (appliesFoldApart(FaceType::Dirichlet)) {
        applyFoldApart(dataTable[dirichletKey], face);
      }
    }
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

void Local::evaluateBatchedTimeDependentBc(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& dataTable,
    SEISSOL_GPU_PARAM recording::ConditionalIndicesTable& indicesTable,
    SEISSOL_GPU_PARAM LTS::Layer& layer,
    SEISSOL_GPU_PARAM double time,
    SEISSOL_GPU_PARAM double timeStepWidth,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {

#ifdef ACL_DEVICE
  using namespace seissol::recording;

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    const ConditionalKey analyticalKey(
        *KernelNames::BoundaryConditions, *ComputationKind::Analytical, face);
    if (indicesTable.find(analyticalKey) != indicesTable.end()) {
      const auto& cellIds =
          indicesTable[analyticalKey].get(inner_keys::Indices::Id::Cells)->getHostData();
      const size_t numElements = cellIds.size();
      auto* analytical =
          reinterpret_cast<real(*)[tensor::INodal::size()]>(layer.var<LTS::AnalyticScratch>());

      runtime.enqueueLoop(
          numElements,
          [this, face, time, timeStepWidth, analytical, &cellIds, &layer](std::size_t index) {
            auto cellId = cellIds.at(index);
            auto data = layer.cellRef(cellId);

            alignas(Alignment) real dofsFaceBoundaryNodal[tensor::INodal::size()];

            assert(initConds_ != nullptr);
            const ApplyAnalyticalSolution applyAnalyticalSolution(initConds_, data);

            analyticalBoundary_.evaluate(data.get<LTS::BoundaryMapping>()[face],
                                         applyAnalyticalSolution,
                                         dofsFaceBoundaryNodal,
                                         time,
                                         timeStepWidth);

            std::memcpy(analytical[index], dofsFaceBoundaryNodal, sizeof(dofsFaceBoundaryNodal));
          });

      auto nodalLfKrnl = deviceNodalLfKrnlPrototype_;
      nodalLfKrnl.INodal = const_cast<const real**>(
          dataTable[analyticalKey].get(inner_keys::Wp::Id::Analytical)->getDeviceDataPtr());
      nodalLfKrnl.AminusT =
          const_cast<const real**>(dataTable[analyticalKey]
                                       .get(inner_keys::Wp::Id::NeighborIntegrationData)
                                       ->getDeviceDataPtr());
      nodalLfKrnl.extraOffset_AminusT =
          SEISSOL_ARRAY_OFFSET(NeighboringIntegrationData, nAmNm1, face);
      nodalLfKrnl.Q =
          dataTable[analyticalKey].get(inner_keys::Wp::Id::FluxTarget)->getDeviceDataPtr();
      nodalLfKrnl.streamPtr = runtime.stream();
      nodalLfKrnl.numElements = numElements;
      nodalLfKrnl.execute(face);
    }
  }
#else
  logError() << "No GPU implementation provided";
#endif // ACL_DEVICE
}

PerformanceEstimate Local::metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes) const {
  PerformanceEstimate estimate;
  estimate += PerformanceEstimate::fromKernel<seissol::kernel::volume>();

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // the trace on every face; both sides of a face cost the same
    estimate += PerformanceEstimate::fromKernel<seissol::kernel::faceTrace>(face, 0);

    if (appliesFoldApart(faceTypes[face])) {
      estimate += PerformanceEstimate::fromKernel<seissol::kernel::foldedBoundaryFlux>(face);
    }

    // Take boundary condition flops into account.
    // Note that this only includes the flops of the kernels but not of the
    // boundary condition implementation.
    // The (probably incorrect) assumption is that they are negligible.
    switch (faceTypes[face]) {
    case FaceType::FreeSurfaceGravity:
      estimate += PerformanceEstimate::fromKernel<seissol::kernel::fsgFlux>(face);
      break;
    case FaceType::Dirichlet:
      estimate += PerformanceEstimate::fromKernel<seissol::kernel::dirichletFlux>(face);
      break;
    case FaceType::Analytical:
      estimate += PerformanceEstimate::fromKernel<seissol::kernel::localFluxNodal>(face);
      estimate +=
          PerformanceEstimate::fromKernel<seissol::kernel::updateINodal>() * ConvergenceOrder;
      break;
    default:
      break;
    }
  }

  // legacy memory estimate
  std::uint64_t reals = 0;

  // star matrices load
  reals += yateto::computeFamilySize<tensor::star>();

  // DOFs write, traces write
  reals += tensor::Q::size() + Cell::NumFaces * tensor::faceTrace::size();

  estimate.bytes = reals * sizeof(real);

  return estimate;
}

} // namespace seissol::kernels::solver::linearck
