/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/

#include "TimeCommon.h"
#include <DataTypes/ConditionalKey.hpp>
#include <DataTypes/ConditionalTable.hpp>
#include <DataTypes/EncodedConstants.hpp>
#include <Initializer/BasicTypedefs.hpp>
#include <Kernels/Time.h>
#include <Kernels/precision.hpp>
#include <Parallel/Runtime/Stream.hpp>
#include <cassert>
#include <stdint.h>
#include <tensor.h>

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
  for (int dofeighbor = 0; dofeighbor < 4; dofeighbor++) {
    assert(((uintptr_t)timeDofs[dofeighbor]) % Alignment == 0);
    assert(((uintptr_t)integrationBuffer[dofeighbor]) % Alignment == 0);
  }
#endif

  /*
   * set/compute time integrated DOFs.
   */
  for (int dofeighbor = 0; dofeighbor < 4; ++dofeighbor) {
    // collect information only in the case that neighboring element contributions are required
    if (faceTypes[dofeighbor] != FaceType::outflow &&
        faceTypes[dofeighbor] != FaceType::dynamicRupture) {
      // check if the time integration is already done (-> copy pointer)
      if ((ltsSetup >> dofeighbor) % 2 == 0) {
        timeIntegrated[dofeighbor] = timeDofs[dofeighbor];
      }
      // integrate the DOFs in time via the derivatives and set pointer to local buffer
      else {
        time.computeIntegral(currentTime[dofeighbor + 1],
                             currentTime[0],
                             currentTime[0] + timeStepWidth,
                             timeDofs[dofeighbor],
                             integrationBuffer[dofeighbor]);

        timeIntegrated[dofeighbor] = integrationBuffer[dofeighbor];
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
    if ((ltsSetup >> (face + 4)) % 2) {
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
  assert(false && "no implementation provided");
#endif
}

} // namespace seissol::kernels
