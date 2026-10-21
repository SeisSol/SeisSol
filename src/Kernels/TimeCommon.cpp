// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "TimeCommon.h"

#include "Common/Constants.h"
#include "Common/Marker.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/LtsSetup.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Parallel/Runtime/Stream.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <stdint.h>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#endif

#ifndef NDEBUG
#include "Alignment.h"

#include <cstdint>
#endif

namespace seissol::kernels {
void TimeCommon::computeNeighborTraces(Time& time,
                                       const LtsSetup& ltsSetup,
                                       const std::array<FaceType, Cell::NumFaces>& faceTypes,
                                       const real* subtimeCoeffs,
                                       const std::array<real*, Cell::NumFaces>& faceNeighbors,
                                       const real* zeroTrace,
                                       const std::array<real*, Cell::NumFaces>& traceBuffers,
                                       std::array<const real*, Cell::NumFaces>& neighborTraces) {
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    if (faceTypes[face] != FaceType::Regular) {
      neighborTraces[face] = zeroTrace;
    } else if (ltsSetup.neighborBuffer(face) != BufferType::Derivatives) {
      // the neighbor wrote the trace over the time step of the cell
      neighborTraces[face] = faceNeighbors[face];
    } else {
      // the neighbor takes a larger time step (see LtsSetup), and the traces of its time
      // derivatives are evaluated over the part of it the cell takes
      assert(!ltsSetup.neighborGTSRelation(face));
      time.evaluateTrace(subtimeCoeffs, faceNeighbors[face], traceBuffers[face]);
      neighborTraces[face] = traceBuffers[face];
    }
  }
}

void TimeCommon::computeBatchedNeighborTraces(
    SEISSOL_GPU_PARAM Time& time,
    SEISSOL_GPU_PARAM const real* subtimeCoeffs,
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  using namespace seissol::recording;
  // the neighbors that take a larger time step, over the part of it the cell takes
  const ConditionalKey key(*KernelNames::NeighborFlux, *ComputationKind::WithLtsDerivatives);
  if (table.find(key) != table.end()) {
    auto& entry = table[key];
    time.evaluateTraceBatched(
        subtimeCoeffs,
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::FaceBuffer))->getDeviceDataPtr()),
        (entry.get(inner_keys::Wp::Id::FaceTraceNeighbor))->getDeviceDataPtr(),
        (entry.get(inner_keys::Wp::Id::FaceTraceNeighbor))->getSize(),
        runtime);
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

} // namespace seissol::kernels
