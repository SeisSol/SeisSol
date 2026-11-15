// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_SIMULATIONFILES_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_SIMULATIONFILES_H_

#include <cstddef>
#include <string>
#include <vector>

namespace seissol::initializer::parameters {

/// Files given one per simulation of a run, with each of them once: what is read from a file serves
/// all simulations that name it.
struct SimulationFiles {
  // the files, in the order in which the simulations first name them
  std::vector<std::string> files;
  // for every simulation, the position of its file in `files`
  std::vector<std::size_t> fileOfSimulation;

  explicit SimulationFiles(const std::vector<std::string>& fileOfEachSimulation) {
    fileOfSimulation.reserve(fileOfEachSimulation.size());
    for (const auto& file : fileOfEachSimulation) {
      std::size_t position = 0;
      while (position < files.size() && files[position] != file) {
        ++position;
      }
      if (position == files.size()) {
        files.push_back(file);
      }
      fileOfSimulation.push_back(position);
    }
  }
};

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_SIMULATIONFILES_H_
