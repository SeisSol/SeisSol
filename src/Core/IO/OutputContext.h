// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_IO_OUTPUTCONTEXT_H_
#define SEISSOL_SRC_CORE_IO_OUTPUTCONTEXT_H_

#include "Core/IO/AsyncIO.h"
#include "Core/Modules/Modules.h"
#include "Core/Parallel/Pin.h"

#include <string>
#include <utils/env.h>

namespace seissol::io {

/// What the outputs use of the run they belong to: the pinning of its threads, its hooks, its
/// asynchronous IO and its environment, and its ID and the prefix of its files. An output reads the
/// last two when it starts; the run sets them when it creates its output directory.
struct OutputContext {
  const parallel::Pinning& pinning;
  Modules& modules;
  AsyncIO& asyncIO;
  utils::Env& env;
  const std::string& runId;
  const std::string& outputPrefix;
};

} // namespace seissol::io

#endif // SEISSOL_SRC_CORE_IO_OUTPUTCONTEXT_H_
