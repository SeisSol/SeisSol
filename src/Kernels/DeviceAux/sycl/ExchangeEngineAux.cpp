// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/DeviceAux/ExchangeEngineAux.h"

#include <cstdint>
#include <utils/logger.h>

namespace seissol::kernels::device::aux::exchangeengine {

// (the engine needs device streams that wait for host memory, which actor-mpi only uses with CUDA
// and HIP)

void request(Transport* /*sender*/,
             std::uint32_t /*send*/,
             Transport* /*receiver*/,
             std::uint32_t /*receive*/,
             void* /*stream*/) {
  logError() << "The exchange engine is only available with CUDA and HIP.";
}

void run(Transport* /*transports*/,
         std::uint32_t /*count*/,
         const std::uint32_t* /*control*/,
         unsigned /*blocks*/,
         void* /*stream*/) {
  logError() << "The exchange engine is only available with CUDA and HIP.";
}

void waitForStop(void* /*stream*/) {}

} // namespace seissol::kernels::device::aux::exchangeengine
