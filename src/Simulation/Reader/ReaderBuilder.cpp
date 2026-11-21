// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ReaderBuilder.h"

#include "Core/Data/DataReader.h"
#include "Simulation/Reader/EasiReader.h"

#include <memory>
#include <string>
#include <utils/logger.h>
#include <vector>

namespace seissol::data {

namespace {

/// The file of a model path, without the prefix `easi:`.
std::string fileOf(const std::string& path) {
  const std::string prefix = "easi:";
  return path.compare(0, prefix.size(), prefix) == 0 ? path.substr(prefix.size()) : path;
}

} // namespace

std::unique_ptr<DataReader> buildReader(const std::string& path,
                                        const std::vector<std::string>& defaultInArgs) {
  logInfo() << "Loading easi file" << fileOf(path);
  return std::make_unique<EasiReader>(fileOf(path), defaultInArgs);
}

} // namespace seissol::data
