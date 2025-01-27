// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_RESULTWRITER_THREADSPINNINGWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_THREADSPINNINGWRITER_H_

#include "Parallel/Pin.h"
#include <string>
#include <vector>

namespace seissol::writer {
class ThreadsPinningWriter {
  public:
  ThreadsPinningWriter(const std::string& outputDirectory) : outputDirectory(outputDirectory) {}
  void write(const seissol::parallel::Pinning& pinning);

  private:
  std::string outputDirectory;
};
} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_THREADSPINNINGWRITER_H_
