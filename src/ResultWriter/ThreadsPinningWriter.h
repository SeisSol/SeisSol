// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_RESULTWRITER_THREADSPINNINGWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_THREADSPINNINGWRITER_H_

#include "Parallel/Pin.h"

#include <string>
#include <utility>
#include <utils/env.h>
#include <vector>

namespace seissol::writer {
class ThreadsPinningWriter {
  public:
  ThreadsPinningWriter(std::string outputDirectory, std::string runId)
      : outputDirectory_(std::move(outputDirectory)), runId_(std::move(runId)) {}
  void write(const seissol::parallel::Pinning& pinning, utils::Env& env);

  private:
  std::string outputDirectory_;
  std::string runId_;
};
} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_THREADSPINNINGWRITER_H_
