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
#include "Initializer/BoundarySetup.h"
#include "Initializer/Typedefs.h"
#include "Kernels/AnalyticalBoundary.h"
#include "Kernels/Common.h"
#include "Kernels/FaceTraces.h"
#include "Kernels/StarOperands.h"
#include "Kernels/WeightAdjustment.h"
#include "Memory/Descriptor/LTS.h"
#include "Model/OperatorLayout.h"
#include "Monitoring/Metric.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <stdint.h>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#endif

GENERATE_HAS_MEMBER(E)
GENERATE_HAS_MEMBER(extraOffset_E)

namespace seissol::kernels::solver::linearckanelastic {

static_assert(!WeightAdjustable || kernels::size<tensor::fluxBuffer>() == tensor::Qext::size(),
              "the faces of a weight-adjusted cell lift into its flux buffer as into the extended "
              "field");

void Local::setGlobalData(const CompoundGlobalData& global) {
  volumeKernelPrototype_.bindGlobals(*global.onHost);
  faceTraceKernelPrototype_.bindGlobals(*global.onHost);
  // the relaxation reads constants too where it is formed at the samples
  localKernelPrototype_.bindGlobals(*global.onHost);

  fsgFlux_.bindGlobals(*global.onHost);
  dirichletFlux_.bindGlobals(*global.onHost);
  foldedBoundaryFlux_.bindGlobals(*global.onHost);
  nodalLfKrnlPrototype_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceFaceTracesKernelPrototype_.bindGlobals(*global.onDevice);
  deviceLocalKernelPrototype_.bindGlobals(*global.onDevice);

  deviceDirichletFlux_.bindGlobals(*global.onDevice);
  deviceFoldedBoundaryFlux_.bindGlobals(*global.onDevice);
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

  // A weight-adjusted cell lifts its face terms without material into its flux buffer, which
  // starts empty here and which the neighbour integral finishes and hands to the field and the
  // memory variables (WeightAdjustable); the volume term goes there through the extended field
  // either way.
  real* fluxTarget = Qext;
  if (WeightAdjustable && data.get<LTS::CellInformation>().weightAdjusted) {
    fluxTarget = data.get<LTS::FluxBuffer>();
    std::fill_n(fluxTarget, kernels::size<tensor::fluxBuffer>(), static_cast<real>(0));
  }

  kernel::volumeExt volKrnl = volumeKernelPrototype_;
  volKrnl.Qext = Qext;
  volKrnl.I = timeIntegratedDoFs;
  kernels::bindStarOperands(volKrnl, data.get<LTS::LocalIntegration>());

  volKrnl.execute();

  // the trace of the time integral on every face, which the corrector of the cell reads together
  // with the one its neighbor writes there
  kernels::writeFaceTraces(faceTraceKernelPrototype_, timeIntegratedDoFs, data);

  const auto& cellBoundaryMapping = data.get<LTS::BoundaryMapping>();
  const auto& materialData = data.get<LTS::Material>();

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // what the condition of the face folds into the local flux solver, where the corrector does
    // not apply it
    if (appliesFoldApart(data.get<LTS::CellInformation>().faceTypes[face])) {
      auto kernel = foldedBoundaryFlux_;
      kernel.Qext = fluxTarget;
      kernel.I = timeIntegratedDoFs;
      kernel.AplusT = data.get<LTS::LocalIntegration>().nApNm1[face];
      kernel.execute(face);
    }

    switch (data.get<LTS::CellInformation>().faceTypes[face]) {
    case FaceType::FreeSurfaceGravity: {
      auto kernel = fsgFlux_;
      kernel.g2m = -2 * this->gravitationalAcceleration_;

      const real localRho = materialData.local->getDensity();
      kernel.rho = &localRho;
      kernel.averageNormalDisplacement = tmp.nodalAvgDisplacements[face].data();

      kernel.Qext = fluxTarget;
      kernel.AminusT = data.get<LTS::NeighboringIntegration>().nAmNm1[face];

      kernel.execute(face);
      break;
    }
    case FaceType::Dirichlet: {
      auto kernel = dirichletFlux_;
      kernel.dirichletOffset = cellBoundaryMapping[face].dirichletOffset;
      kernel.dt = timeStepWidth;

      kernel.Qext = fluxTarget;
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
      nodalLfKrnl.Qext = fluxTarget;
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

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // the trace on every face; both sides of a face cost the same
    estimate += PerformanceEstimate::fromKernel<seissol::kernel::faceTrace>(face, 0);
    if (appliesFoldApart(faceTypes[face])) {
      estimate += PerformanceEstimate::fromKernel<seissol::kernel::foldedBoundaryFlux>(face);
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

  // legacy memory estimate
  std::uint64_t reals = 0;

  // star matrices load
  reals += yateto::computeFamilySize<tensor::star>() + tensor::w::size() + tensor::W::size() +
           tensor::E::size();

  // DOFs write, traces write
  reals += tensor::Q::size() + tensor::Qane::size() + Cell::NumFaces * tensor::faceTrace::size();

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

  // A weight-adjusted cell lifts its face terms into its flux buffer, which starts empty here and
  // which the neighbour integral finishes (WeightAdjustable); the volume term goes to the field
  // and the memory variables through the extended field either way.
  WeightAdjustmentT::clearBatched(dataTable, runtime);

  // The predictor has written the volume term to the extended field already
  // (Spacetime::computeBatchedAder).
  const ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  if (dataTable.find(key) != dataTable.end()) {
    auto& entry = dataTable[key];
    const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();

    kernel::gpu_local localKrnl = deviceLocalKernelPrototype_;
    auto tmpMem = runtime.memoryHandle<real>(
        (yateto::getMaxTmpMemRequired(localKrnl, deviceDirichletFlux_, deviceFoldedBoundaryFlux_) *
         numElements) /
        sizeof(real));

    // the trace of the time integral on every face, which the corrector of the cell reads together
    // with the one its neighbor writes there
    kernels::writeFaceTracesBatched(deviceFaceTracesKernelPrototype_, dataTable, runtime);

    // The Dirichlet condition of the boundary faces goes where the faces of a cell lift into -- the
    // extended field, which the local kernel below hands on to the field and the memory variables,
    // or the flux buffer of a weight-adjusted cell; the free surface with gravity is not defined
    // for the materials of this solver.
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      const ConditionalKey dirichletKey(
          *KernelNames::BoundaryConditions, *ComputationKind::Dirichlet, face);
      if (dataTable.find(dirichletKey) == dataTable.end()) {
        continue;
      }
      auto& bcEntry = dataTable[dirichletKey];
      const auto bcElements = bcEntry.get(inner_keys::Wp::Id::Dofs)->getSize();
      auto** fluxTargets = bcEntry.get(inner_keys::Wp::Id::FluxTarget)->getDeviceDataPtr();

      auto bcKernel = deviceDirichletFlux_;
      bcKernel.numElements = bcElements;
      bcKernel.dirichletOffset = const_cast<const real**>(
          bcEntry.get(inner_keys::Wp::Id::DirichletOffset)->getDeviceDataPtr());
      bcKernel.dt = timeStepWidth;
      bcKernel.Qext = fluxTargets;
      bcKernel.AminusT = const_cast<const real**>(
          bcEntry.get(inner_keys::Wp::Id::NeighborIntegrationData)->getDeviceDataPtr());
      bcKernel.extraOffset_AminusT = SEISSOL_ARRAY_OFFSET(NeighboringIntegrationData, nAmNm1, face);
      bcKernel.linearAllocator.initialize(tmpMem.get());
      bcKernel.streamPtr = runtime.stream();
      bcKernel.execute(face);

      // what the condition folds into the local flux solver, where the corrector does not apply
      // it
      if constexpr (appliesFoldApart(FaceType::Dirichlet)) {
        auto foldKernel = deviceFoldedBoundaryFlux_;
        foldKernel.numElements = bcElements;
        foldKernel.Qext = fluxTargets;
        foldKernel.I =
            const_cast<const real**>(bcEntry.get(inner_keys::Wp::Id::Idofs)->getDeviceDataPtr());
        foldKernel.AplusT = const_cast<const real**>(
            bcEntry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
        foldKernel.extraOffset_AplusT = SEISSOL_ARRAY_OFFSET(LocalIntegrationData, nApNm1, face);
        foldKernel.linearAllocator.initialize(tmpMem.get());
        foldKernel.streamPtr = runtime.stream();
        foldKernel.execute(face);
      }
    }

    localKrnl.numElements = numElements;
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
    // where the material varies inside the cell, the relaxation is formed
    // from what it says at the sample points and the kernel takes no
    // matrix at all
    set_E(localKrnl,
          const_cast<const real**>(
              entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr()));
    set_extraOffset_E(localKrnl, SEISSOL_OFFSET(LocalIntegrationData, specific.E));
    kernels::bindSourceOperandsBatched(
        localKrnl,
        const_cast<const real**>(
            entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr()));
    localKrnl.linearAllocator.initialize(tmpMem.get());
    localKrnl.streamPtr = runtime.stream();

    SEISSOL_OFFSET_ASSERT(LocalIntegrationData, specific.W);
    SEISSOL_OFFSET_ASSERT(LocalIntegrationData, specific.w);
    SEISSOL_OFFSET_ASSERT(LocalIntegrationData, specific.E);

    localKrnl.execute();
  }

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
