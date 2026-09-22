// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Local.h"

#include "Common/Constants.h"
#include "Common/Marker.h"
#include "Common/Offset.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Interface.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <utils/logger.h>
#include <yateto.h>

namespace seissol::kernels::solver::nonlinearck {

void Local::setGlobalData(const CompoundGlobalData& global) {
  cellIntegral_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceCellIntegral_.bindGlobals(*global.onDevice);
#endif
}

void Local::computeIntegral(real* timeIntegratedDoFs,
                            LTS::Ref& data,
                            LocalTmp& tmp,
                            double /*time*/,
                            double /*timeStepWidth*/) {
  // What the cell owes itself: the volume term and the source of the
  // quantities that carry no flux. Both are constant maps on what the
  // predictor produced.
  kernel::damageCellIntegral krnl = cellIntegral_;
  krnl.I = timeIntegratedDoFs;
  krnl.Q = data.get<LTS::Dofs>();
  krnl.sourceI = tmp.sourceIntegral;
  krnl.materialParameters = data.get<LTS::LocalIntegration>().specific.parameters;
  // The reference direction a derivative is taken along is not a physical one:
  // which physical direction it is, is the geometry the cell was meshed with,
  // and the volume term needs it as much as the recursion does.
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    krnl.star(i) = data.get<LTS::LocalIntegration>().starMatrices[i];
  }
  krnl.execute();

  // Every face flux is assembled where both wave speeds are, which is the
  // neighbouring integration -- including the faces that have no neighbour,
  // whose ghost rule is folded into their own pair of matrices at setup.
}

void Local::computeBatchedIntegral(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& dataTable,
    SEISSOL_GPU_PARAM recording::ConditionalMaterialTable& materialTable,
    SEISSOL_GPU_PARAM recording::ConditionalIndicesTable& indicesTable,
    SEISSOL_GPU_PARAM double timeStepWidth,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  using namespace seissol::recording;

  const ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  if (dataTable.find(key) == dataTable.end()) {
    return;
  }
  auto& entry = dataTable[key];

  kernel::gpu_damageCellIntegral krnl = deviceCellIntegral_;
  krnl.numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
  krnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
  krnl.I = const_cast<const real**>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
  krnl.sourceI = const_cast<const real**>(
      (entry.get(inner_keys::Wp::Id::SourceIntegrals))->getDeviceDataPtr());

  constexpr auto ParametersOffset =
      offsetof(LocalIntegrationData, specific) + offsetof(NonLinearLocalData, parameters);
  static_assert(ParametersOffset % sizeof(real) == 0,
                "The material of a cell is not aligned to the real size.");
  const auto** localIntegrationPtrs = const_cast<const real**>(
      (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
  krnl.materialParameters = localIntegrationPtrs;
  krnl.extraOffset_materialParameters = ParametersOffset / sizeof(real);

  // As on the host: the geometry of the cell, reached the way the material is.
  SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData, starMatrices);
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    krnl.star(i) = localIntegrationPtrs;
    krnl.extraOffset_star(i) = SEISSOL_ARRAY_OFFSET(LocalIntegrationData, starMatrices, i);
  }

  auto tmpMem = runtime.memoryHandle<real>((yateto::getMaxTmpMemRequired(krnl) * krnl.numElements) /
                                           sizeof(real));
  krnl.linearAllocator.initialize(tmpMem.get());
  krnl.streamPtr = runtime.stream();
  krnl.execute();
#else
  logError() << "No GPU implementation provided";
#endif
}

void Local::evaluateBatchedTimeDependentBc(
    recording::ConditionalPointersToRealsTable& /*dataTable*/,
    recording::ConditionalIndicesTable& /*indicesTable*/,
    LTS::Layer& /*layer*/,
    double /*time*/,
    double /*timeStepWidth*/,
    seissol::parallel::runtime::StreamRuntime& /*runtime*/) {
  // Nothing to do, as on the host: this material evaluates no time-dependent
  // boundary condition. Its local integral is the volume and the source term,
  // which computeBatchedIntegral has already added, and every face -- a
  // boundary face too -- is assembled in the neighboring integration. The
  // launch code calls this right after the cell integral, so running that
  // integral here again added the volume and the source term twice per step.
}

PerformanceEstimate
    Local::metrics(const std::array<FaceType, Cell::NumFaces>& /*faceTypes*/) const {
  auto estimate = PerformanceEstimate::fromKernel<kernel::damageCellIntegral>();

  std::uint64_t reals = 0;
  // the transported tensor and the source integral in, the state in and out
  reals += tensor::I::size() + tensor::sourceI::size() + 2 * tensor::Q::size();

  estimate.bytes = reals * sizeof(real);
  return estimate;
}

} // namespace seissol::kernels::solver::nonlinearck
