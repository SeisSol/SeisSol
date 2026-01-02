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
void TimeCommon::computeIntegrals(Time& time,
                                  const LtsSetup& ltsSetup,
                                  const std::array<FaceType, Cell::NumFaces>& faceTypes,
                                  const real* timeCoeffs,
                                  const real* subtimeCoeffs,
                                  real* const timeDofs[4],
                                  real integrationBuffer[4][tensor::I::size()],
                                  real* timeIntegrated[4]) {
  // call the more general assembly
  /*
   * assert valid input.
   */

#ifndef NDEBUG
  // alignment of the time derivatives/integrated dofs and the buffer
  for (std::size_t dofneighbor = 0; dofneighbor < Cell::NumFaces; dofneighbor++) {
    assert(reinterpret_cast<uintptr_t>(timeDofs[dofneighbor]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(integrationBuffer[dofneighbor]) % Alignment == 0);
  }
#endif

  /*
   * set/compute time integrated DOFs.
   */
  for (std::size_t dofneighbor = 0; dofneighbor < Cell::NumFaces; ++dofneighbor) {
    // collect information only in the case that neighboring element contributions are required
    if (faceTypes[dofneighbor] != FaceType::Outflow &&
        faceTypes[dofneighbor] != FaceType::DynamicRupture) {
      // check if the time integration is already done (-> copy pointer)
      if (!ltsSetup.neighborHasDerivatives(dofneighbor)) {
        timeIntegrated[dofneighbor] = timeDofs[dofneighbor];
      }
      // integrate the DOFs in time via the derivatives and set pointer to local buffer
      else {
        // select the time coefficients next: dependent on if we have GTS as a neighbor, use GTS
        // coefficients (timeCoeffs); otherwise use LTS coefficients (subtimeCoeffs) for a large
        // time cluster neighbor. IMPORTANT: make sure to not land in this code path if the neighbor
        // time cluster is less than the local time cluster. It shouldn't happen with the current
        // setup; but just be aware of it when changing things. In that case, enforce the "GTS
        // relation" instead; then everything will work again.

        const auto* coeffs = ltsSetup.neighborGTSRelation(dofneighbor) ? timeCoeffs : subtimeCoeffs;
        time.evaluate(coeffs, timeDofs[dofneighbor], integrationBuffer[dofneighbor]);

        timeIntegrated[dofneighbor] = integrationBuffer[dofneighbor];
      }
    }
  }
}

void TimeCommon::computeBatchedIntegrals(
    SEISSOL_GPU_PARAM Time& time,
    SEISSOL_GPU_PARAM const real* timeCoeffs,
    SEISSOL_GPU_PARAM const real* subtimeCoeffs,
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  using namespace seissol::recording;
  // Compute time integrated dofs using neighbors derivatives using the GTS relation,
  // i.e. the expansion point is around 'timeStepStart'
  ConditionalKey key(*KernelNames::NeighborFlux, *ComputationKind::WithGtsDerivatives);
  if (table.find(key) != table.end()) {
    auto& entry = table[key];
    time.evaluateBatched(
        timeCoeffs,
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr()),
        (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr(),
        (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
        runtime);
  }

  // Compute time integrated dofs using neighbors derivatives using the LTS relation,
  // i.e. the expansion point is around '0'
  key = ConditionalKey(*KernelNames::NeighborFlux, *ComputationKind::WithLtsDerivatives);
  if (table.find(key) != table.end()) {
    auto& entry = table[key];
    time.evaluateBatched(
        subtimeCoeffs,
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr()),
        (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr(),
        (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
        runtime);
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

} // namespace seissol::kernels
