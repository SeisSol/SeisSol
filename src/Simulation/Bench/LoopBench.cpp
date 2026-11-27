// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "LoopBench.h"

#include "Simulation/Bench/Proxy.h"

#include <memory>

namespace seissol::proxy {

auto makeLoopBench(Mode mode, const BenchSetup& setup) -> std::unique_ptr<LoopBench> {
  return mode == Mode::Cluster ? makeClusterBench(setup) : makeKernelBench(setup);
}

} // namespace seissol::proxy
