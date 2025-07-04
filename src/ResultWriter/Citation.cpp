// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Citation.h"
#include <Parallel/MPI.h>
#include <fstream>
#include <string>
#include <unordered_map>

namespace {
const std::unordered_map<seissol::CitationFeature, std::string> CitationData{
    {seissol::CitationFeature::Core, R"(
TODO
)"}};
} // namespace

namespace seissol {

void CitationHandler::enable(seissol::CitationFeature feature) { features.push_back(feature); }
void CitationHandler::print(const std::string& prefix) {
  if (MPI::mpi.rank() == 0) {
    std::ofstream ofile(prefix + "-citation.bib");
    for (const auto& feature : features) {
      ofile << CitationData.at(feature);
    }
  }
}

} // namespace seissol
