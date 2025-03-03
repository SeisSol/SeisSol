// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_RESULTWRITER_MINISEISSOLWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_MINISEISSOLWRITER_H_

#include <string>
#include <vector>

namespace seissol::writer {
class MiniSeisSolWriter {
  public:
  MiniSeisSolWriter(const char* outputDirectory) : outputDirectory(outputDirectory) {}
  void write(double elapsedTime, double weight);

  private:
  std::string outputDirectory;
};
} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_MINISEISSOLWRITER_H_
