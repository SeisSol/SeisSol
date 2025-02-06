// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "TimeCommon.h"
#include <DataTypes/ConditionalTable.h>
#include <Initializer/BasicTypedefs.h>
#include <Kernels/Precision.h>
#include <Kernels/Time.h>
#include <Parallel/Runtime/Stream.h>
#include <cassert>
#include <stdint.h>
#include <tensor.h>

#include "utils/logger.h"

#ifdef ACL_DEVICE
#include <DataTypes/ConditionalKey.h>
#include <DataTypes/EncodedConstants.h>
#endif

#ifndef NDEBUG
#include "Common/Constants.h"
#include <cstdint>
#endif

namespace seissol::kernels {

void TimeCommon::computeIntegrals(Time& time,
                                  unsigned short ltsSetup,
                                  const FaceType faceTypes[4],
                                  const double currentTime[5],
                                  double timeStepWidth,
                                  real* const timeDofs[4],
                                  real integrationBuffer[4][tensor::I::size()],
                                  real* timeIntegrated[4]) {
  /*
   * assert valid input.
   */
  // only lower 10 bits are used for lts encoding
  assert(ltsSetup < 2048);

#ifndef NDEBUG
  // alignment of the time derivatives/integrated dofs and the buffer
  for (int dofneighbor = 0; dofneighbor < 4; dofneighbor++) {
    assert(reinterpret_cast<uintptr_t>(timeDofs[dofneighbor]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(integrationBuffer[dofneighbor]) % Alignment == 0);
  }
#endif

  /*
   * set/compute time integrated DOFs.
   */
  for (int dofneighbor = 0; dofneighbor < 4; ++dofneighbor) {
    // collect information only in the case that neighboring element contributions are required
    if (faceTypes[dofneighbor] != FaceType::Outflow &&
        faceTypes[dofneighbor] != FaceType::DynamicRupture) {
      // check if the time integration is already done (-> copy pointer)
      if ((ltsSetup >> dofneighbor) % 2 == 0) {
        timeIntegrated[dofneighbor] = timeDofs[dofneighbor];
      }
      // integrate the DOFs in time via the derivatives and set pointer to local buffer
      else {
        time.computeIntegral(currentTime[dofneighbor + 1],
                             currentTime[0],
                             currentTime[0] + timeStepWidth,
                             timeDofs[dofneighbor],
                             integrationBuffer[dofneighbor]);

        timeIntegrated[dofneighbor] = integrationBuffer[dofneighbor];
      }
    }
  }
}

void TimeCommon::computeIntegrals(Time& time,
                                  unsigned short ltsSetup,
                                  const FaceType faceTypes[4],
                                  const double timeStepStart,
                                  const double timeStepWidth,
                                  real* const timeDofs[4],
                                  real integrationBuffer[4][tensor::I::size()],
                                  real* timeIntegrated[4]) {
  double startTimes[5];
  startTimes[0] = timeStepStart;
  startTimes[1] = startTimes[2] = startTimes[3] = startTimes[4] = 0;

  // adjust start times for GTS on derivatives
  for (unsigned int face = 0; face < 4; face++) {
    if (((ltsSetup >> (face + 4)) % 2) != 0) {
      startTimes[face + 1] = timeStepStart;
    }
  }

  // call the more general assembly
  computeIntegrals(time,
                   ltsSetup,
                   faceTypes,
                   startTimes,
                   timeStepWidth,
                   timeDofs,
                   integrationBuffer,
                   timeIntegrated);
}

void TimeCommon::computeBatchedIntegrals(Time& time,
                                         const double timeStepStart,
                                         const double timeStepWidth,
                                         ConditionalPointersToRealsTable& table,
                                         seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  // Compute time integrated dofs using neighbours derivatives using the GTS relation,
  // i.e. the expansion point is around 'timeStepStart'
  ConditionalKey key(*KernelNames::NeighborFlux, *ComputationKind::WithGtsDerivatives);
  if (table.find(key) != table.end()) {
    auto& entry = table[key];
    time.computeBatchedIntegral(
        timeStepStart,
        timeStepStart,
        timeStepStart + timeStepWidth,
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr()),
        (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr(),
        (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
        runtime);
  }

  // Compute time integrated dofs using neighbours derivatives using the LTS relation,
  // i.e. the expansion point is around '0'
  key = ConditionalKey(*KernelNames::NeighborFlux, *ComputationKind::WithLtsDerivatives);
  if (table.find(key) != table.end()) {
    auto& entry = table[key];
    time.computeBatchedIntegral(
        0.0,
        timeStepStart,
        timeStepStart + timeStepWidth,
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
