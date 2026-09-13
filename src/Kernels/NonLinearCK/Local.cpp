// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Local.h"

#include "Common/Constants.h"
#include "Common/Marker.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Interface.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

#include <array>
#include <cstdint>

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
  krnl.rhoInv = data.get<LTS::LocalIntegration>().specific.parameters.rhoInv;
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
  logError() << "No GPU implementation provided";
}

void Local::evaluateBatchedTimeDependentBc(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& dataTable,
    SEISSOL_GPU_PARAM recording::ConditionalIndicesTable& indicesTable,
    SEISSOL_GPU_PARAM LTS::Layer& layer,
    SEISSOL_GPU_PARAM double time,
    SEISSOL_GPU_PARAM double timeStepWidth,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
  logError() << "No GPU implementation provided";
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
